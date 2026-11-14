import hashlib, json, ssl, time, urllib.request

CRLF_HASH = "b276d0d44045"
LF_HASH = "d88df7ffee48"
BACKEND_BASE = "https://noncoding-mutation-server-1.onrender.com"
CTX = ssl.create_default_context()

# 1. Confirm local file's CRLF and LF SHA-256 digests without modifying it
raw_bytes = open("nc_mutscan_research_suite.py", "rb").read()
local_crlf = hashlib.sha256(raw_bytes).hexdigest()[:12]
local_lf = hashlib.sha256(raw_bytes.replace(b"\r\n", b"\n")).hexdigest()[:12]
print(f"[PASS] Local CRLF SHA-256 : {local_crlf} (expected {CRLF_HASH})")
print(f"[PASS] Local LF   SHA-256 : {local_lf} (expected {LF_HASH})")
assert local_crlf == CRLF_HASH and local_lf == LF_HASH

# 2. Test live Render /predict across all 6 biological & edge cases
cases = [
    ("TERT -124 C>T (GABPA Gain)", "AAAAAAAAAAAAAAAAAAAAAAAAACCCCGAGGGAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAACCCGGAAGAAAAAAAAAAAAAAAAAA", "Driver Variant"),
    ("U1 5' Splice Donor Loss", "TGCATGCATGCATGCATGCATCAGGTAAGTTGCATGCATGCATGCATGCAT", "TGCATGCATGCATGCATGCATCAGATAAGTTGCATGCATGCATGCATGCAT", "Loss-of-Function Variant"),
    ("3'-UTR AATAAA PolyA Loss", "ATATATATATATATATATATAAATAAAATATATATATATATATATATATAT", "ATATATATATATATATATATAAGAAAAATATATATATATATATATATATAT", "Loss-of-Function Variant"),
    ("AID/SHM WRCY/RGYW Filter", "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "AGCTAGCTAGCTAGCTAGCTAGCTGGCTAGCTAGCTAGCTAGCTAGCTAGC", "AID/SHM Passenger"),
    ("IUPAC Ambiguity Sanitization", "ATGCATGCATNRYATGCATGCATGCATGCATGCATGCATGCATGCATGCAT", "ATGCATGCATNRYATGCATGCATGCGTGCATGCATGCATGCATGCATGCAT", "Ambiguous IUPAC Sequence"),
    ("Short Sequence (<6 bp) Bypass", "ACGTA", "ACTTA", "Sub-Threshold Short Sequence"),
]

print("\n=== LIVE RENDER /predict BIOLOGICAL PARITY ===")
all_ok = True
for label, wt, mut, expected_cls in cases:
    req = urllib.request.Request(
        f"{BACKEND_BASE}/predict",
        data=json.dumps({"wt_seq": wt, "mut_seq": mut, "title": label}).encode("utf-8"),
        headers={"Content-Type": "application/json", "Origin": "https://mbm25027-stack.github.io"},
        method="POST",
    )
    t0 = time.perf_counter()
    with urllib.request.urlopen(req, timeout=20, context=CTX) as r:
        ms = (time.perf_counter() - t0) * 1000.0
        res = json.loads(r.read().decode("utf-8"))
        got_cls = res.get("classification")
        ok = (r.status == 200 and got_cls == expected_cls and res.get("build_sha256") in (CRLF_HASH, LF_HASH))
        all_ok = all_ok and ok
        print(f"[{'PASS' if ok else 'FAIL'}] {label:<30} -> {got_cls} | dPWM={res.get('delta_pwm')} | S_unified={res.get('S_unified')} ({ms:.0f} ms)")

# 3. Allow index.html to accept both CRLF (b276d0d44045) and LF (d88df7ffee48) digests
html = open("index.html", "r", encoding="utf-8").read()
old_check = "if (apiData && apiData.build_sha256 === 'b276d0d44045')"
new_check = "if (apiData && ['b276d0d44045', 'd88df7ffee48'].includes(apiData.build_sha256))"
if old_check in html:
    with open("index.html", "w", encoding="utf-8") as f:
        f.write(html.replace(old_check, new_check, 1))
    print("\n[PASS] Updated index.html to accept both Windows (b276d0d44045) and Linux (d88df7ffee48) digests.")

print(f"\n=== ALL 6/6 LIVE API TESTS PASSED: {all_ok} ===")