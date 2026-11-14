#!/usr/bin/env python3
"""
server_benchmark.py - End-to-end latency, cold-start, concurrency, determinism and robustness
checks for the NC-MutScan HTTP API. Standard library only.

The request body is a JSON template with {placeholders} filled from --variants (TSV: chrom pos ref alt
[sequence]) or reused as-is. Example:
  python server_benchmark.py --url https://noncoding-mutation-server-1.onrender.com/predict \
     --template '{"chrom":"{chrom}","pos":{pos},"ref":"{ref}","alt":"{alt}"}' \
     --variants variants.tsv --n 200 --concurrency 1 4 16 --out report.json
Adjust --template to the server's real schema. Add --expect-field build_hash to check the response
carries the frozen-build identifier. Use only against a server you own.
"""
import argparse, concurrent.futures as cf, hashlib, json, statistics, sys, time
import urllib.request, urllib.error

def call(url, body, method="POST", headers=None, timeout=60):
    data = body.encode() if isinstance(body, str) else body
    req = urllib.request.Request(url, data=data, method=method,
                                 headers=headers if headers is not None else {"Content-Type": "application/json"})
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read(), time.perf_counter() - t0
    except urllib.error.HTTPError as e:
        return e.code, e.read(), time.perf_counter() - t0
    except Exception as e:  # timeout / connection reset
        return -1, str(e).encode(), time.perf_counter() - t0

def pct(xs, p):
    xs = sorted(xs); k = (len(xs) - 1) * p / 100
    lo, hi = int(k), min(int(k) + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)
def summarize(lat):
    ms = [x * 1000 for x in lat]
    return {"n": len(ms), "mean_ms": round(statistics.mean(ms), 2), "p50_ms": round(pct(ms, 50), 2),
            "p95_ms": round(pct(ms, 95), 2), "p99_ms": round(pct(ms, 99), 2), "max_ms": round(max(ms), 2)}

def load_variants(path):
    out = []
    for line in open(path):
        if line.strip() and not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            out.append(dict(zip(["chrom", "pos", "ref", "alt", "sequence"], f)))
    return out
def make_bodies(template, variants, n):
    if not variants: return [template] * n
    def fill(v):
        s = template
        for k in ("chrom", "pos", "ref", "alt", "sequence"):
            s = s.replace("{" + k + "}", v.get(k) or "")
        return s
    return [fill(variants[i % len(variants)]) for i in range(n)]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--url", required=True); ap.add_argument("--template", required=True)
    ap.add_argument("--variants"); ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--concurrency", type=int, nargs="+", default=[1, 4, 16])
    ap.add_argument("--timeout", type=float, default=90)
    ap.add_argument("--expect-field", default="", help="JSON field that must be present (e.g. build_hash)")
    ap.add_argument("--out", default="server_report.json")
    a = ap.parse_args()
    variants = load_variants(a.variants) if a.variants else []
    bodies = make_bodies(a.template, variants, a.n)
    rep = {"url": a.url, "time": time.strftime("%Y-%m-%d %H:%M:%S")}

    # 1) cold start: first request (run this after the server has been idle, e.g. >15 min)
    st, body, dt = call(a.url, bodies[0], timeout=a.timeout)
    rep["cold_start"] = {"status": st, "seconds": round(dt, 2), "likely_cold": dt > 3.0}
    print(f"[cold] status={st} {dt:.2f}s" + ("  <-- likely cold start" if dt > 3 else ""))
    if st != 200: print("First request failed; check --template/--url.", body[:200]); sys.exit(1)
    try:
        j = json.loads(body); rep["response_fields"] = sorted(j) if isinstance(j, dict) else "non-object"
        if a.expect_field: rep["has_expected_field"] = a.expect_field in j; print(f"[schema] {a.expect_field}: {rep['has_expected_field']}")
    except ValueError: rep["response_fields"] = "non-JSON"

    # 2) warm sequential latency
    lat, errs = [], 0
    for b in bodies:
        s, _, d = call(a.url, b, timeout=a.timeout)
        if s == 200: lat.append(d)
        else: errs += 1
    rep["warm_sequential"] = {**summarize(lat), "errors": errs} if lat else {"errors": errs}
    print("[warm]", rep["warm_sequential"])

    # 3) concurrency sweep
    rep["concurrency"] = {}
    for c in a.concurrency:
        t0 = time.perf_counter(); lats = []; bad = 0
        with cf.ThreadPoolExecutor(max_workers=c) as ex:
            for s, _, d in ex.map(lambda b: call(a.url, b, timeout=a.timeout), bodies):
                (lats.append(d) if s == 200 else None); bad += s != 200
        wall = time.perf_counter() - t0
        rep["concurrency"][c] = {**(summarize(lats) if lats else {}), "errors": bad,
                                 "throughput_rps": round(len(bodies) / wall, 1)}
        print(f"[conc={c}]", rep["concurrency"][c])

    # 4) determinism: identical request -> identical response bytes
    hashes = {hashlib.sha256(call(a.url, bodies[0], timeout=a.timeout)[1]).hexdigest() for _ in range(20)}
    rep["deterministic"] = len(hashes) == 1; print("[determinism] identical across 20 calls:", rep["deterministic"])

    # 5) robustness: bad inputs must yield 4xx (never 5xx / hang) and the server must stay alive
    big = json.dumps({"sequence": "A" * 2_000_000})
    cases = {"empty_body": ("", "POST", None), "invalid_json": ("{not json", "POST", None),
             "empty_object": ("{}", "POST", None), "null_json": ("null", "POST", None),
             "oversize_2MB": (big, "POST", None), "null_bytes": (bodies[0].replace("A", "A\\u0000", 1), "POST", None),
             "unicode": (bodies[0].replace("A", "Å", 1), "POST", None),
             "wrong_content_type": (bodies[0], "POST", {"Content-Type": "text/plain"}),
             "GET_on_predict": ("", "GET", {}), "sql_like": (bodies[0].replace("A", "A';--", 1), "POST", None)}
    rep["robustness"] = {}
    for name, (b, m, h) in cases.items():
        s, rb, d = call(a.url, b if m == "POST" else None, method=m, headers=h, timeout=a.timeout)
        ok = s != -1 and s < 500
        rep["robustness"][name] = {"status": s, "ok_no_5xx": ok, "seconds": round(d, 2)}
        print(f"[robust] {name:20s} status={s:<4} {'OK' if ok else 'FAIL'}")
    s, _, _ = call(a.url, bodies[0], timeout=a.timeout)
    rep["alive_after_abuse"] = s == 200; print("[alive after abuse]", rep["alive_after_abuse"])

    json.dump(rep, open(a.out, "w"), indent=2); print("Report written:", a.out)

if __name__ == "__main__":
    main()