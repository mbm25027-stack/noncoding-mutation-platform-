import json, os, ssl, subprocess, time, urllib.error, urllib.request

FRONTEND_URL = "https://mbm25027-stack.github.io/noncoding-mutation-server/"
BACKEND_BASE = "https://noncoding-mutation-server-1.onrender.com"
PREDICT_URL = f"{BACKEND_BASE}/predict"
CTX = ssl.create_default_context()

def call(url, method="GET", payload=None, headers=None, timeout=30):
    h = {"User-Agent": "NC-MutScan-Verifier/1.1", "Accept": "application/json, text/html, */*"}
    if headers:
        h.update(headers)
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    if data is not None:
        h["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=h, method=method)
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
            return r.status, dict(r.headers), r.read().decode("utf-8", "replace"), (time.perf_counter() - t0) * 1000
    except urllib.error.HTTPError as e:
        return e.code, dict(e.headers), e.read().decode("utf-8", "replace"), (time.perf_counter() - t0) * 1000
    except Exception as e:
        return 0, {}, str(e), (time.perf_counter() - t0) * 1000

print("=== 1. LOCAL GIT & INDEX.HTML CHECK ===")
for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if d not in {"venv", ".venv", ".git", "__pycache__"}]
    for f in files:
        if f.lower() in {"index.html", "app.py", "main.py", "server.py"}:
            print("Found:", os.path.join(root, f))
try:
    rem = subprocess.check_output(["git", "remote", "-v"], stderr=subprocess.STDOUT, text=True).strip()
    print("Git remotes:\n" + rem)
except Exception as e:
    print("Git remote check:", e)

print("\n=== 2. FRONTEND CHECK ===")
st, _, body, ms = call(FRONTEND_URL)
print(f"GET {FRONTEND_URL} -> HTTP {st} ({ms:.0f} ms)")

print("\n=== 3. POLLING RENDER UNTIL CONTAINER BOOTS (MAX 90s) ===")
awake = False
for attempt in range(1, 19):
    st, hdrs, body, ms = call(f"{BACKEND_BASE}/health", timeout=25)
    if st == 404:
        st, hdrs, body, ms = call(BACKEND_BASE, timeout=25)
    loading = "Application loading" in body
    print(f"Attempt {attempt:02d}: HTTP {st} ({ms:.0f} ms) | loading_screen={loading}")
    if st in (200, 404, 405) and not loading:
        awake = True
        break
    time.sleep(5)

if not awake:
    print("[FAIL] Render did not finish booting within 90s. Check Render dashboard Logs for startup error.")
    raise SystemExit(1)

_, chdrs, _, _ = call(PREDICT_URL, method="OPTIONS", headers={
    "Origin": "https://mbm25027-stack.github.io",
    "Access-Control-Request-Method": "POST",
    "Access-Control-Request-Headers": "Content-Type",
})
print("CORS Allow-Origin:", chdrs.get("Access-Control-Allow-Origin") or chdrs.get("access-control-allow-origin", "MISSING"))

print("\n=== 4. SCHEMA & BIOLOGICAL ACCURACY CHECKS ===")
wt_tert = "CCCCTCCCGGGTCCCCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCC"
mut_tert = "CCCCTCCCGGGTCCCCGGCCCAGCCCCTTCCGGGCCCTCCCAGCCCCTCCC"
schemas = [
    ("ref_seq/mut_seq", lambda w, m: {"ref_seq": w, "mut_seq": m}),
    ("wt_seq/alt_seq", lambda w, m: {"wt_seq": w, "alt_seq": m}),
    ("wild_type/mutant", lambda w, m: {"wild_type": w, "mutant": m}),
    ("sequence/mutated_sequence", lambda w, m: {"sequence": w, "mutated_sequence": m}),
    ("ref/alt", lambda w, m: {"ref": w, "alt": m}),
    ("wt/mut", lambda w, m: {"wt": w, "mut": m}),
]
builder = None
for name, fn in schemas:
    st, _, body, ms = call(PREDICT_URL, method="POST", payload=fn(wt_tert, mut_tert))
    if st == 200:
        builder = fn
        print(f"[PASS] Schema matched: {name} ({ms:.0f} ms)")
        print("Raw TERT -124C>T response:", body[:350])
        break
    else:
        print(f"Tried {name} -> HTTP {st}: {body[:140]}")

if builder:
    cases = [
        ("TC-01 TERT -124C>T", wt_tert, mut_tert),
        ("TC-06 TP53 Sp1 loss", "TTTTTGGGCGGGGCCGGGCGGGGCCGGGCGGGGCCGGGCGGGGCCAAAAAA", "TTTTTGGGCGGGGCCGGGCGGGGCCGGGCGAGGCCGGGCGGGGCCAAAAAA"),
        ("TC-03/12 AID WRCY/RGYW", "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC", "AGCTAGCTAGCTAGCTAGCTAGCTGGCTAGCTAGCTAGCTAGCTAGCTAGC"),
        ("U1 5' Splice Donor Loss", "TGCATGCATGCATGCATGCATCAGGTAAGTTGCATGCATGCATGCATGCAT", "TGCATGCATGCATGCATGCATCAGATAAGTTGCATGCATGCATGCATGCAT"),
        ("3'-UTR AATAAA PolyA Loss", "ATATATATATATATATATATAAATAAAATATATATATATATATATATATAT", "ATATATATATATATATATATAAGAAAAATATATATATATATATATATATAT"),
        ("TC-04 IUPAC N/R/Y", "ATGCATGCATNRYATGCATGCATGCATGCATGCATGCATGCATGCATGCAT", "ATGCATGCATNRYATGCATGCATGCGTGCATGCATGCATGCATGCATGCAT"),
        ("TC-10 Short <6bp", "ACGTA", "ACTTA"),
    ]
    for label, w, m in cases:
        st, _, body, ms = call(PREDICT_URL, method="POST", payload=builder(w, m))
        print(f"[{st}] {label:<26} ({ms:.0f} ms) -> {body[:200]}")