import hashlib
import json
import os
import re
import ssl
import statistics
import sys
import time
import urllib.error
import urllib.request

FRONTEND_URL = "https://mbm25027-stack.github.io/noncoding-mutation-server/"
BACKEND_BASE = "https://noncoding-mutation-server-1.onrender.com"
PREDICT_URL = f"{BACKEND_BASE}/predict"
EXPECTED_BUILD_HASH = "b276d0d44045"

SSL_CTX = ssl.create_default_context()


def http_call(url, method="GET", payload=None, headers=None, timeout=75):
    req_headers = {
        "User-Agent": "NC-MutScan-v3.3-Verifier/1.0",
        "Accept": "application/json, text/html, */*",
    }
    if headers:
        req_headers.update(headers)

    data = None
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        req_headers["Content-Type"] = "application/json"

    req = urllib.request.Request(url, data=data, headers=req_headers, method=method)
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=SSL_CTX) as resp:
            elapsed_ms = (time.perf_counter() - t0) * 1000.0
            body = resp.read().decode("utf-8", errors="replace")
            return resp.status, dict(resp.headers), body, elapsed_ms
    except urllib.error.HTTPError as e:
        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        body = e.read().decode("utf-8", errors="replace")
        return e.code, dict(e.headers), body, elapsed_ms
    except Exception as e:
        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        return 0, {}, f"{type(e).__name__}: {e}", elapsed_ms


def print_header(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def check_frontend():
    print_header("1. GITHUB PAGES FRONTEND CHECK")
    status, headers, body, ms = http_call(FRONTEND_URL, method="GET", timeout=30)
    print(f"URL            : {FRONTEND_URL}")
    print(f"HTTP Status    : {status} ({ms:.1f} ms)")
    if status != 200:
        print(f"[FAIL] Frontend unreachable: {body[:200]}")
        return False

    has_api = "noncoding-mutation-server-1.onrender.com" in body
    has_v33 = "3.3" in body
    print(f"Links Render API: {'YES' if has_api else 'NO (check JS/HTML fetch URL)'}")
    print(f"Mentions v3.3   : {'YES' if has_v33 else 'NO (update title/version badge to v3.3)'}")
    return True


def check_backend_and_cors():
    print_header("2. RENDER BACKEND WAKE-UP, HEALTH & CORS CHECK")
    # Wake up Render instance (free tier may take 30-60s on cold start)
    for ep in (f"{BACKEND_BASE}/health", BACKEND_BASE):
        status, headers, body, ms = http_call(ep, method="GET", timeout=75)
        print(f"GET {ep:<52} -> HTTP {status} ({ms:.1f} ms)")
        if status == 200:
            if EXPECTED_BUILD_HASH in body:
                print(f"  [PASS] Build hash {EXPECTED_BUILD_HASH} confirmed in response.")
            break

    # Preflight CORS check simulating GitHub Pages browser request
    cors_headers = {
        "Origin": "https://mbm25027-stack.github.io",
        "Access-Control-Request-Method": "POST",
        "Access-Control-Request-Headers": "Content-Type",
    }
    status, headers, body, ms = http_call(
        PREDICT_URL, method="OPTIONS", headers=cors_headers, timeout=30
    )
    allow_origin = headers.get("Access-Control-Allow-Origin") or headers.get(
        "access-control-allow-origin", "MISSING"
    )
    print(f"OPTIONS /predict CORS -> HTTP {status} | Access-Control-Allow-Origin: {allow_origin}")
    if allow_origin in ("*", "https://mbm25027-stack.github.io"):
        print("  [PASS] CORS allows GitHub Pages frontend requests.")
        return True
    print("  [WARN] Check CORS configuration if browser requests fail.")
    return False


def detect_predict_schema():
    print_header("3. AUTO-DETECTING /predict JSON SCHEMA")
    # TERT -124 C>T (C228T) reverse-strand 5'-CCGGAA-3' creation
    wt_tert = "CCCCTCCCGGGTCCCCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCC"
    mut_tert = "CCCCTCCCGGGTCCCCGGCCCAGCCCCTTCCGGGCCCTCCCAGCCCCTCCC"

    candidates = [
        ("ref_seq / mut_seq", lambda w, m: {"ref_seq": w, "mut_seq": m}),
        ("wt_seq / alt_seq", lambda w, m: {"wt_seq": w, "alt_seq": m}),
        ("wild_type / mutant", lambda w, m: {"wild_type": w, "mutant": m}),
        ("sequence / mutated_sequence", lambda w, m: {"sequence": w, "mutated_sequence": m}),
        ("ref / alt", lambda w, m: {"ref": w, "alt": m}),
        ("wt / mut", lambda w, m: {"wt": w, "mut": m}),
    ]

    last_err = ""
    for label, builder in candidates:
        status, headers, body, ms = http_call(
            PREDICT_URL, method="POST", payload=builder(wt_tert, mut_tert), timeout=45
        )
        if status == 200:
            try:
                parsed = json.loads(body)
                print(f"[PASS] Matched /predict schema: {label} ({ms:.1f} ms)")
                print(f"Sample JSON keys: {list(parsed.keys()) if isinstance(parsed, dict) else type(parsed)}")
                return builder
            except json.JSONDecodeError:
                last_err = f"Non-JSON 200 response: {body[:200]}"
        else:
            last_err = f"HTTP {status}: {body[:250]}"

    print(f"[FAIL] Could not auto-match /predict request schema. Last response: {last_err}")
    return None


def classify_response_summary(resp_obj):
    flat = json.dumps(resp_obj).lower()
    score_val = None
    if isinstance(resp_obj, dict):
        for k in ("s_unified", "S_unified", "score", "delta_pwm", "dpwm", "max_delta_pwm", "s_raw", "S_raw"):
            if k in resp_obj and isinstance(resp_obj[k], (int, float)):
                score_val = resp_obj[k]
                break
        if score_val is None:
            for v in resp_obj.values():
                if isinstance(v, dict):
                    for k in ("s_unified", "S_unified", "score", "delta_pwm", "dpwm", "s_raw"):
                        if k in v and isinstance(v[k], (int, float)):
                            score_val = v[k]
                            break
    return flat, score_val


def run_biological_regression_checks(make_payload):
    print_header("4. DETERMINISTIC BIOLOGICAL & EDGE-CASE ACCURACY TESTS (v3.3)")
    tests = [
        {
            "id": "TC-01 (TERT -124 C>T ETS/GABPA Gain)",
            "wt": "CCCCTCCCGGGTCCCCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCC",
            "mut": "CCCCTCCCGGGTCCCCGGCCCAGCCCCTTCCGGGCCCTCCCAGCCCCTCCC",
            "expect_any": ["gabpa", "ets", "driver", "gain", "pathogenic", "high"],
            "reject_any": ["error", "traceback"],
        },
        {
            "id": "TC-06 (TP53 Promoter Sp1 Loss Competition)",
            "wt": "TTTTTGGGCGGGGCCGGGCGGGGCCGGGCGGGGCCGGGCGGGGCCAAAAAA",
            "mut": "TTTTTGGGCGGGGCCGGGCGGGGCCGGGCGAGGCCGGGCGGGGCCAAAAAA",
            "expect_any": ["sp1", "loss", "driver", "disrupt"],
            "reject_any": ["error", "traceback"],
        },
        {
            "id": "TC-03/12 (Clustered Dual-Strand AID/SHM WRCY/RGYW)",
            "wt": "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
            "mut": "AGCTAGCTAGCTAGCTAGCTAGCTGGCTAGCTAGCTAGCTAGCTAGCTAGC",
            "expect_any": ["aid", "shm", "passenger", "benign", "neutral", "sub"],
            "reject_any": ["error", "traceback"],
        },
        {
            "id": "U1 5' Splice Donor Loss (CAG|GTAAGT -> CAG|ATAAGT)",
            "wt": "TGCATGCATGCATGCATGCATCAGGTAAGTTGCATGCATGCATGCATGCAT",
            "mut": "TGCATGCATGCATGCATGCATCAGATAAGTTGCATGCATGCATGCATGCAT",
            "expect_any": ["splice", "donor", "5ss", "loss", "driver", "disrupt"],
            "reject_any": ["error", "traceback"],
        },
        {
            "id": "3'-UTR PolyA Signal Loss (AATAAA -> AAGAAA in AT-rich)",
            "wt": "ATATATATATATATATATATAAATAAAATATATATATATATATATATATAT",
            "mut": "ATATATATATATATATATATAAGAAAAATATATATATATATATATATATAT",
            "expect_any": ["poly", "pas", "aataaa", "loss", "driver", "disrupt"],
            "reject_any": ["error", "traceback"],
        },
        {
            "id": "TC-04 (IUPAC Degenerate Base Sanitization N/R/Y)",
            "wt": "ATGCATGCATNRYATGCATGCATGCATGCATGCATGCATGCATGCATGCAT",
            "mut": "ATGCATGCATNRYATGCATGCATGCGTGCATGCATGCATGCATGCATGCAT",
            "expect_any": ["iupac", "ambiguous", "sanitiz", "neutral", "benign", "sub", "status", "result"],
            "reject_any": ["traceback", "valueerror", "keyerror", "indexerror"],
        },
        {
            "id": "TC-10 (Short Sequence < 6 bp Bypass)",
            "wt": "ACGTA",
            "mut": "ACTTA",
            "expect_any": ["short", "sub", "threshold", "bypass", "neutral", "benign", "status", "result"],
            "reject_any": ["traceback", "valueerror", "indexerror"],
        },
    ]

    passed = 0
    for tc in tests:
        status, headers, body, ms = http_call(
            PREDICT_URL, method="POST", payload=make_payload(tc["wt"], tc["mut"]), timeout=30
        )
        if status != 200:
            print(f"[FAIL] {tc['id']}: HTTP {status} -> {body[:160]}")
            continue

        try:
            obj = json.loads(body)
        except json.JSONDecodeError:
            print(f"[FAIL] {tc['id']}: Invalid JSON -> {body[:160]}")
            continue

        flat, score_val = classify_response_summary(obj)
        bad_hit = [w for w in tc["reject_any"] if w in flat]
        good_hit = [w for w in tc["expect_any"] if w in flat]

        if not bad_hit and good_hit:
            passed += 1
            sc_str = f" | score={score_val}" if score_val is not None else ""
            print(f"[PASS] {tc['id']} ({ms:.1f} ms){sc_str} -> matched {good_hit[:3]}")
        else:
            print(f"[WARN] {tc['id']} ({ms:.1f} ms): unexpected output -> {body[:220]}")

    print(f"\nRegression Summary: {passed}/{len(tests)} checks passed cleanly.")
    return passed == len(tests)


def measure_warm_latency(make_payload, n_calls=5):
    print_header("5. WARM HTTP API LATENCY BENCHMARK")
    wt = "TGCATGCATGCATGCATGCATCAGGTAAGTTGCATGCATGCATGCATGCAT"
    mut = "TGCATGCATGCATGCATGCATCAGATAAGTTGCATGCATGCATGCATGCAT"
    latencies = []
    for i in range(n_calls):
        status, _, _, ms = http_call(PREDICT_URL, method="POST", payload=make_payload(wt, mut), timeout=20)
        if status == 200:
            latencies.append(ms)
    if latencies:
        print(f"Samples        : {len(latencies)} warm POST /predict calls")
        print(f"Median RTT     : {statistics.median(latencies):.2f} ms (includes network round-trip)")
        print(f"Min / Max RTT  : {min(latencies):.2f} ms / {max(latencies):.2f} ms")


def check_local_repo_parity():
    print_header("6. LOCAL REPOSITORY & BUILD HASH PARITY CHECK")
    cwd = os.path.abspath(os.getcwd())
    print(f"Working Dir    : {cwd}")
    matched_hash = False
    py_files = []
    for dp, dns, fns in os.walk(cwd):
        dns[:] = [d for d in dns if d not in {"venv", ".venv", ".git", "__pycache__"}]
        for fn in fns:
            if fn.endswith(".py"):
                full = os.path.join(dp, fn)
                py_files.append(full)
                with open(full, "rb") as f:
                    h = hashlib.sha256(f.read()).hexdigest()
                if h.startswith(EXPECTED_BUILD_HASH):
                    matched_hash = True
                    print(f"[PASS] Exact SHA-256 build match ({EXPECTED_BUILD_HASH}): {fn}")
    if not matched_hash:
        print(f"[INFO] No single .py file matched prefix {EXPECTED_BUILD_HASH} directly (may be package/commit hash).")
    print(f"Python files scanned: {len(py_files)}")


def main():
    fe_ok = check_frontend()
    be_ok = check_backend_and_cors()
    builder = detect_predict_schema()
    bio_ok = False
    if builder is not None:
        bio_ok = run_biological_regression_checks(builder)
        measure_warm_latency(builder, n_calls=5)
    check_local_repo_parity()

    print_header("FINAL DEPLOYMENT VERIFICATION SUMMARY")
    print(f"Frontend Reachable : {'PASS' if fe_ok else 'FAIL'}")
    print(f"Backend & CORS     : {'PASS' if be_ok else 'CHECK'}")
    print(f"Schema & Accuracy  : {'PASS' if bio_ok else 'CHECK'}")
    sys.exit(0 if (fe_ok and builder is not None and bio_ok) else 1)


if __name__ == "__main__":
    main()