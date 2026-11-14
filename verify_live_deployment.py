import hashlib, http.client, json, ssl, time, urllib.error, urllib.request

EXPECTED_HASH = "b276d0d44045"
FRONTEND_URL = "https://mbm25027-stack.github.io/noncoding-mutation-platform-/"
BACKEND_BASE = "https://noncoding-mutation-server-1.onrender.com"
CTX = ssl.create_default_context()

def call(url, method="GET", payload=None, headers=None, timeout=25):
    h = {
        "User-Agent": "NC-MutScan-Live-Verifier/3.3",
        "Accept": "application/json, text/html, */*",
    }
    if headers:
        h.update(headers)
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    if data is not None:
        h["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=h, method=method)
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
            return r.status, dict(r.headers), r.read().decode("utf-8", "replace"), (time.perf_counter() - t0) * 1000.0
    except urllib.error.HTTPError as e:
        try:
            err_body = e.read().decode("utf-8", "replace")
        except http.client.IncompleteRead as ir:
            err_body = (ir.partial or b"").decode("utf-8", "replace")
        except Exception:
            err_body = ""
        return e.code, dict(e.headers), err_body, (time.perf_counter() - t0) * 1000.0
    except http.client.IncompleteRead as ir:
        return 503, {}, (ir.partial or b"").decode("utf-8", "replace"), (time.perf_counter() - t0) * 1000.0
    except Exception as e:
        return 0, {}, str(e), (time.perf_counter() - t0) * 1000.0

print("=== 1. LIVE GITHUB PAGES FRONTEND CHECK ===")
st, _, html, ms = call(FRONTEND_URL)
has_v33 = "NC-MutScan v3.3" in html
has_sha = EXPECTED_HASH in html
has_api = "noncoding-mutation-server-1.onrender.com" in html
print(f"GET {FRONTEND_URL} -> HTTP {st} ({ms:.0f} ms)")
print(f"[{'PASS' if (st == 200 and has_v33 and has_sha and has_api) else 'CHECK'}] v3.3_label={has_v33} | sha256_tag={has_sha} | render_hook={has_api}")

print("\n=== 2. LIVE RENDER BACKEND HEALTH & COLD-START CHECK (MAX 120s) ===")
health_ok = False
for attempt in range(1, 25):
    st, hdrs, body, ms = call(f"{BACKEND_BASE}/health")
    if st == 200 and "Application loading" not in body:
        try:
            data = json.loads(body)
            ver = data.get("version")
            sha = data.get("build_sha256")
            print(f"[{'PASS' if sha == EXPECTED_HASH else 'CHECK'}] GET /health -> HTTP 200 ({ms:.0f} ms) | version={ver} | build_sha256={sha}")
            health_ok = (sha == EXPECTED_HASH)
            break
        except Exception:
            pass
    print(f"Attempt {attempt:02d}: HTTP {st} ({ms:.0f} ms) - container booting/deploying...")
    time.sleep(5)

if not health_ok:
    print("[INFO] Render did not return build_sha256=b276d0d44045 yet.")
    print("       Verify in dashboard.render.com -> noncoding-mutation-server-1 -> Settings:")
    print("       Repository = mbm25027-stack/noncoding-mutation-platform- (branch: main, startCommand: python server.py)")
    raise SystemExit(1)

print("\n=== 3. LIVE CORS & /predict BIOLOGICAL PARITY CHECKS ===")
_, chdrs, _, _ = call(f"{BACKEND_BASE}/predict", method="OPTIONS", headers={
    "Origin": "https://mbm25027-stack.github.io",
    "Access-Control-Request-Method": "POST",
    "Access-Control-Request-Headers": "Content-Type",
})
cors_origin = chdrs.get("Access-Control-Allow-Origin") or chdrs.get("access-control-allow-origin", "MISSING")
print(f"[{'PASS' if cors_origin in ('*', 'https://mbm25027-stack.github.io') else 'CHECK'}] CORS Allow-Origin: {cors_origin}")

live_cases = [
    ("TERT -124 C>T (GABPA Gain)", "AAAAAAAAAAAAAAAAAAAAAAAAACCCCGAGGGAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAACCCGGAAGAAAAAAAAAAAAAAAAAA", "Driver Variant"),
    ("U1 5' Splice Donor Loss", "TGCATGCATGCATGCATGCATCAGGTAAGTTGCATGCATGCATGCATGCAT", "TGCATGCATGCATGCATGCATCAGATAAGTTGCATGCATGCATGCATGCAT", "Loss-of-Function Variant"),
    ("3'-UTR AATAAA PolyA Loss", "ATATATATATATATATATATAAATAAAATATATATATATATATATATATAT", "ATATATATATATATATATATAAGAAAAATATATATATATATATATATATAT", "Loss-of-Function Variant"),
    ("AID/SHM WRCY/RGYW Filter", "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "AGCTAGCTAGCTAGCTAGCTAGCTGGCTAGCTAGCTAGCTAGCTAGCTAGC", "AID/SHM Passenger"),
    ("IUPAC Ambiguity Sanitization", "ATGCATGCATNRYATGCATGCATGCATGCATGCATGCATGCATGCATGCAT", "ATGCATGCATNRYATGCATGCATGCGTGCATGCATGCATGCATGCATGCAT", "Ambiguous IUPAC Sequence"),
    ("Short Sequence (<6 bp) Bypass", "ACGTA", "ACTTA", "Sub-Threshold Short Sequence"),
]

all_passed = True
for label, wt, mut, expected_cls in live_cases:
    st, _, body, ms = call(f"{BACKEND_BASE}/predict", method="POST", payload={"wt_seq": wt, "mut_seq": mut, "title": label})
    res = json.loads(body) if st == 200 else {}
    got_cls = res.get("classification")
    ok = (st == 200 and got_cls == expected_cls and res.get("build_sha256") == EXPECTED_HASH)
    all_passed = all_passed and ok
    print(f"[{'PASS' if ok else 'FAIL'}] {label:<30} -> {got_cls} | dPWM={res.get('delta_pwm')} | S_unified={res.get('S_unified')} ({ms:.0f} ms)")

local_sha = hashlib.sha256(open("nc_mutscan_research_suite.py", "rb").read()).hexdigest()[:12]
print(f"\n=== FINAL LIVE DEPLOYMENT STATUS ===")
print(f"Local Frozen Core SHA-256 : {local_sha}")
print(f"Frontend + Backend Parity : {'VERIFIED' if (st == 200 and has_v33 and has_sha and health_ok and all_passed) else 'CHECK OUTPUT ABOVE'}")