import hashlib, io, json

EXPECTED_HASH = "b276d0d44045"
srv = open("server.py", "r", encoding="utf-8").read()

wsgi_block = '''
def _wsgi_response(start_response, status_line, body_bytes, content_type="application/json"):
    headers = [
        ("Content-Type", content_type),
        ("Content-Length", str(len(body_bytes))),
        ("Access-Control-Allow-Origin", "*"),
        ("Access-Control-Allow-Methods", "GET, POST, OPTIONS"),
        ("Access-Control-Allow-Headers", "Content-Type, Authorization"),
    ]
    start_response(status_line, headers)
    return [body_bytes]


def app(environ, start_response):
    """PEP 3333 WSGI entrypoint for Gunicorn (gunicorn server:app)."""
    method = environ.get("REQUEST_METHOD", "GET").upper()
    path = environ.get("PATH_INFO", "/")

    if method == "OPTIONS":
        return _wsgi_response(start_response, "200 OK", b"")

    if method == "GET":
        if path in ("/", "/health", "/api/health"):
            resp = {
                "status": "online",
                "engine": "NC-MutScan v3.3",
                "version": "NC-MutScan v3.3",
                "build_sha256": CORE_SHORT_SHA256,
                "full_sha256": CORE_FULL_SHA256,
                "build_verified": CORE_SHORT_SHA256 == EXPECTED_BUILD_PREFIX,
                "endpoints": ["/health", "/predict", "/api/analyze"],
            }
            return _wsgi_response(start_response, "200 OK", json.dumps(resp, indent=2).encode("utf-8"))

        if path == "/index.html":
            index_file = os.path.join(BASE_DIR, "index.html")
            if os.path.exists(index_file):
                with open(index_file, "rb") as f:
                    return _wsgi_response(start_response, "200 OK", f.read(), content_type="text/html; charset=utf-8")

        return _wsgi_response(start_response, "404 Not Found", json.dumps({"error": f"Endpoint not found: {path}"}).encode("utf-8"))

    if method == "POST":
        try:
            raw_len = str(environ.get("CONTENT_LENGTH") or "0").strip()
            content_length = int(raw_len) if raw_len.isdigit() else 0
            if content_length > 2 * 1024 * 1024:
                return _wsgi_response(start_response, "413 Payload Too Large", json.dumps({"error": "Payload exceeds 2 MB limit."}).encode("utf-8"))

            body_bytes = environ["wsgi.input"].read(content_length) if content_length > 0 else b""
            post_data = body_bytes.decode("utf-8") if body_bytes else "{}"
            payload = json.loads(post_data) if post_data.strip() else {}

            if path in ("/api/analyze_sv", "/analyze_sv"):
                if sv_engine and hasattr(sv_engine, "analyze_structural_variant"):
                    res = sv_engine.analyze_structural_variant(
                        payload.get("sv_type", "TRANSLOCATION"),
                        payload.get("chrom_a", ""),
                        int(payload.get("pos_a", 0)),
                        payload.get("chrom_b", ""),
                        int(payload.get("pos_b", 0)),
                    )
                else:
                    res = {"error": "Structural Variant engine unavailable."}
                return _wsgi_response(start_response, "200 OK", json.dumps(res, indent=2).encode("utf-8"))

            if path in ("/", "/predict", "/api/predict", "/api/analyze", "/analyze"):
                result = run_v33_prediction(payload)
                if result is None:
                    err = {
                        "error": "Provide 'wt_seq' and 'mut_seq' (or 'fasta' with header >chr:pos | Name (REF>ALT)).",
                        "version": "NC-MutScan v3.3",
                        "build_sha256": CORE_SHORT_SHA256,
                    }
                    return _wsgi_response(start_response, "400 Bad Request", json.dumps(err).encode("utf-8"))
                return _wsgi_response(start_response, "200 OK", json.dumps(result, indent=2).encode("utf-8"))

            return _wsgi_response(start_response, "404 Not Found", json.dumps({"error": f"Endpoint not found: {path}"}).encode("utf-8"))
        except json.JSONDecodeError:
            return _wsgi_response(start_response, "400 Bad Request", json.dumps({"error": "Invalid JSON payload."}).encode("utf-8"))
        except Exception as e:
            return _wsgi_response(start_response, "500 Internal Server Error", json.dumps({"error": str(e)}).encode("utf-8"))

    return _wsgi_response(start_response, "405 Method Not Allowed", json.dumps({"error": f"Method not allowed: {method}"}).encode("utf-8"))
'''

if "def app(environ, start_response):" not in srv:
    srv = srv.replace("def run():", wsgi_block + "\n\ndef run():", 1)
    with open("server.py", "w", encoding="utf-8") as f:
        f.write(srv)
    print("[PASS] Added PEP 3333 WSGI 'app' callable to server.py.")

# Verify WSGI app callable directly
import server
assert callable(getattr(server, "app", None)), "server.app missing!"

captured = {}
def fake_start_response(status, headers):
    captured["status"] = status
    captured["headers"] = dict(headers)

# 1. Test GET /health via WSGI
out = b"".join(server.app({"REQUEST_METHOD": "GET", "PATH_INFO": "/health"}, fake_start_response))
health = json.loads(out.decode("utf-8"))
assert captured["status"] == "200 OK" and health["build_sha256"] == EXPECTED_HASH
print(f"[PASS] WSGI GET /health  -> {captured['status']} | version={health['version']} | sha256={health['build_sha256']}")

# 2. Test POST /predict via WSGI
body = json.dumps({
    "wt_seq": "AAAAAAAAAAAAAAAAAAAAAAAAACCCCGAGGGAAAAAAAAAAAAAAAAA",
    "mut_seq": "AAAAAAAAAAAAAAAAAAAAAAAAACCCGGAAGAAAAAAAAAAAAAAAAAA",
    "title": "TERT -124 C>T",
}).encode("utf-8")
out = b"".join(server.app({
    "REQUEST_METHOD": "POST",
    "PATH_INFO": "/predict",
    "CONTENT_LENGTH": str(len(body)),
    "wsgi.input": io.BytesIO(body),
}, fake_start_response))
pred = json.loads(out.decode("utf-8"))
assert captured["status"] == "200 OK" and pred["classification"] == "Driver Variant" and pred["build_sha256"] == EXPECTED_HASH
print(f"[PASS] WSGI POST /predict -> {captured['status']} | {pred['classification']} | dPWM={pred['delta_pwm']} | S_unified={pred['S_unified']}")

core_sha = hashlib.sha256(open("nc_mutscan_research_suite.py", "rb").read()).hexdigest()[:12]
assert core_sha == EXPECTED_HASH
print(f"[VERIFIED] Frozen core engine untouched: sha256:{core_sha}")