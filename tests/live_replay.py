#!/usr/bin/env python3
"""
live_replay.py -- replay the golden regression set (dev) and the blind subset against the LIVE API.
Standard library only. Exit code 1 on any failure.

  python tests/live_replay.py --base-url https://your-host --golden tests/golden/golden_v3.3.0.json --release RELEASE.json

ASSUMPTION TO CHECK: the prediction route is POST <base>/predict with JSON {"wt_seq","mut_seq","flank_bp"}
and the reply carries S_raw and classification (top level or under "result"). Change --predict-path or
dig() if your API differs (see verify_live_parity.py).
"""
import argparse, json, sys, time, urllib.error, urllib.request

NEG = ("benign", "neutral", "polymorphism", "passenger", "sub-threshold", "subthreshold", "non-driver", "compet")
POS = ("driver", "pathogenic", "loss", "gain", "splice", "polya", "disrupt")


def frozen_call(cls):   # same rule as is_frozen_driver_call in the frozen scripts
    t = str(cls or "").lower()
    if any(n in t for n in NEG):
        return 0
    return 1 if any(p in t for p in POS) else 0


def http(method, url, body=None, timeout=60, retries=3):
    data = None if body is None else json.dumps(body).encode()
    req = urllib.request.Request(url, data=data, method=method,
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            return e.code, e.read()
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise SystemExit(f"Cannot reach {url}: {last}")


def dig(res, key):
    if isinstance(res, dict):
        if key in res:
            return res[key]
        for v in res.values():
            if isinstance(v, dict) and key in v:
                return v[key]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-url", required=True)
    ap.add_argument("--golden", required=True)
    ap.add_argument("--release", default="RELEASE.json")
    ap.add_argument("--predict-path", default="/predict")
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--sleep", type=float, default=0.0)
    ap.add_argument("--out", default="parity_report.json")
    a = ap.parse_args()
    base = a.base_url.rstrip("/")
    rel = json.load(open(a.release, encoding="utf-8"))
    gold = json.load(open(a.golden, encoding="utf-8"))
    fails, report = [], {"base_url": base, "checks": {}, "sets": {}}

    def check(name, ok, detail=""):
        report["checks"][name] = {"ok": bool(ok), "detail": detail}
        print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
        if not ok:
            fails.append(name)

    st, _ = http("GET", base + "/health")
    check("health 200", st == 200, f"(status {st})")

    st, body = http("GET", base + "/version")
    try:
        v = json.loads(body) if st == 200 else {}
    except ValueError:
        v = {}
    check("/version 200", st == 200 and bool(v), f"(status {st})")
    check("version.git_tag matches RELEASE.json", v.get("git_tag") == rel["git_tag"], f"(live {v.get('git_tag')})")
    check("version.core hash matches pin", rel["core_sha256"].startswith(str(v.get("core_short_sha256") or "~")),
          f"(live {v.get('core_short_sha256')})")
    check("version.build_verified true", v.get("build_verified") is True)
    check("golden file made on this core", gold["meta"]["core_sha256"] == rel["core_sha256"])

    for name in ("dev", "blind"):
        recs = [r for r in gold["records"] if r["set"] == name]
        mism, bad_score, bad_call, errs = [], 0, 0, 0
        live_conf = {"tp": 0, "fp": 0, "tn": 0, "fn": 0}
        for r in recs:
            st, body = http("POST", base + a.predict_path,
                            {"wt_seq": r["wt_seq"], "mut_seq": r["mut_seq"], "flank_bp": gold["meta"]["flank_bp"]})
            try:
                res = json.loads(body)
            except ValueError:
                res = None
            s, cls = dig(res, "S_raw"), dig(res, "classification")
            if st != 200 or s is None:
                errs += 1
                mism.append((r["id"], f"http {st}, no S_raw"))
                continue
            call = frozen_call(cls)
            if abs(float(s) - r["S_raw"]) > a.tol:
                bad_score += 1
                mism.append((r["id"], f"S_raw live {s} vs golden {r['S_raw']}"))
            if call != r["call"]:
                bad_call += 1
                mism.append((r["id"], f"call live {call} vs golden {r['call']}"))
            if r.get("label") in (0, 1):
                key = ("tp" if call else "fn") if r["label"] == 1 else ("fp" if call else "tn")
                live_conf[key] += 1
            if a.sleep:
                time.sleep(a.sleep)
        ok = bool(recs) and not mism
        check(f"replay {name} (n={len(recs)})", ok,
              f"(score mismatches {bad_score}, call mismatches {bad_call}, errors {errs})")
        for m in mism[:10]:
            print("   ", m)
        report["sets"][name] = {"n": len(recs), "score_mismatch": bad_score, "call_mismatch": bad_call,
                                "errors": errs, "live_confusion": live_conf}

    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"\n{'ALL CHECKS PASSED' if not fails else 'FAILED: ' + ', '.join(fails)}  (report: {a.out})")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
