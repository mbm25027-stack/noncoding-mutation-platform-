import os, re, ssl, urllib.error, urllib.request

CTX = ssl.create_default_context()
ALT_FRONTEND = "https://mbm25027-stack.github.io/noncoding-mutation-platform-/"
BACKEND = "https://noncoding-mutation-server-1.onrender.com"

def get_url(url):
    req = urllib.request.Request(url, headers={"User-Agent": "NC-MutScan-Inspector/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=20, context=CTX) as r:
            return r.status, r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode("utf-8", "replace")
    except Exception as e:
        return 0, str(e)

print("=== 1. CHECK ACTUAL GITHUB PAGES URL ===")
st, body = get_url(ALT_FRONTEND)
print(f"GET {ALT_FRONTEND} -> HTTP {st} (len={len(body)})")
if st == 200:
    print("Title snippet:", re.findall(r"<title>(.*?)</title>", body, re.I | re.S))

print("\n=== 2. INSPECT LOCAL SERVER & FRONTEND ROUTES ===")
for fname in ("app.py", "main.py", "server.py", "Procfile", "render.yaml", "index.html", r"files (1)\index.html"):
    if os.path.isfile(fname):
        txt = open(fname, "r", encoding="utf-8", errors="ignore").read()
        print(f"\n--- {fname} ({len(txt)} chars) ---")
        for line in txt.splitlines():
            if any(k in line for k in ("@app.", "@router.", "route(", "fetch(", "onrender.com", "nc_mutscan", "gunicorn", "uvicorn", "Mutated_Donors")):
                print("  ", line.strip()[:120])

print("\n=== 3. PROBE LIVE RENDER ENDPOINTS ===")
for path in ("/openapi.json", "/docs", "/scan", "/scan_sequence", "/v3", "/api/scan", "/analyze"):
    st, body = get_url(BACKEND + path)
    print(f"GET {path:<18} -> HTTP {st} | {body[:140].replace(chr(10), ' ')}")