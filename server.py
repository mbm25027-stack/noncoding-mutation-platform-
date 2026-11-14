import hashlib
import http.server
import importlib.util
import json
import os
import re
import sys
from urllib.parse import urlparse

PORT = int(os.environ.get("PORT", 8000))
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CORE_PATH = os.path.join(BASE_DIR, "nc_mutscan_research_suite.py")
EXPECTED_BUILD_PREFIX = "b276d0d44045"


def load_frozen_v33_engine():
    with open(CORE_PATH, "rb") as f:
        full_sha256 = hashlib.sha256(f.read()).hexdigest()
    spec = importlib.util.spec_from_file_location("nc_mutscan_v33", CORE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod, full_sha256


NC_ENGINE, CORE_FULL_SHA256 = load_frozen_v33_engine()
CORE_SHORT_SHA256 = CORE_FULL_SHA256[:12]

# Optional legacy DB / SV support (kept intact if present)
NONCODING_DB = {}
MODEL_PATH = os.path.join(BASE_DIR, "trained_server_model", "noncoding_driver_db.json")
if os.path.exists(MODEL_PATH):
    try:
        with open(MODEL_PATH, "r", encoding="utf-8") as f:
            NONCODING_DB = json.load(f)
    except Exception:
        NONCODING_DB = {}

try:
    import sv_engine
except Exception:
    sv_engine = None


def parse_single_fasta_to_pair(fasta_text):
    lines = [ln.strip() for ln in str(fasta_text or "").strip().splitlines() if ln.strip()]
    header = ""
    seq_parts = []
    for ln in lines:
        if ln.startswith(">"):
            header = ln[1:].strip()
        else:
            seq_parts.append(re.sub(r"\s+", "", ln).upper())
    mut_seq = "".join(seq_parts)
    parts = [p.strip() for p in header.split("|")]
    lm = re.search(r"(chr[\dXYM]+:\d+)", header, re.I)
    locus = lm.group(1) if lm else (parts[0] if parts and parts[0] else "Custom Locus")
    title = parts[1].replace("_", " ") if len(parts) >= 2 and parts[1] else (header or "Variant Analysis")

    vm = re.search(r"([ACGT])\s*>\s*([ACGT])", header, re.I)
    wt_seq = None
    flank_bp = 25
    if vm and len(mut_seq) >= 6:
        ref_b = vm.group(1).upper()
        alt_b = vm.group(2).upper()
        mid = len(mut_seq) // 2
        # Search outward from midpoint for alt_b to reconstruct exact WT sequence
        best_idx = None
        for offset in range(len(mut_seq)):
            for cand in (mid - offset, mid + offset):
                if 0 <= cand < len(mut_seq) and mut_seq[cand] == alt_b:
                    best_idx = cand
                    break
            if best_idx is not None:
                break
        if best_idx is not None:
            wt_seq = mut_seq[:best_idx] + ref_b + mut_seq[best_idx + 1:]
            flank_bp = best_idx

    return wt_seq, mut_seq, flank_bp, locus, title


def run_v33_prediction(payload):
    wt_seq = (
        payload.get("wt_seq")
        or payload.get("ref_seq")
        or payload.get("wild_type")
        or payload.get("wt")
        or payload.get("ref")
        or ""
    ).strip()
    mut_seq = (
        payload.get("mut_seq")
        or payload.get("alt_seq")
        or payload.get("mutant")
        or payload.get("mut")
        or payload.get("alt")
        or payload.get("mutated_sequence")
        or ""
    ).strip()

    raw_flank = payload.get("flank_bp", None)
    flank_bp = int(raw_flank) if raw_flank is not None else 25
    locus = str(payload.get("locus", "Custom Locus")).strip() or "Custom Locus"
    title = str(payload.get("title", "Variant Analysis")).strip() or "Variant Analysis"

    # Support single FASTA string or 'sequence' field from web UI
    fasta_raw = (payload.get("fasta") or payload.get("sequence") or "").strip()
    if (not wt_seq or not mut_seq) and fasta_raw:
        parsed_wt, parsed_mut, parsed_flank, parsed_locus, parsed_title = parse_single_fasta_to_pair(fasta_raw)
        if parsed_wt and parsed_mut:
            wt_seq = parsed_wt
            mut_seq = parsed_mut
            flank_bp = int(payload.get("flank_bp", parsed_flank))
            if locus == "Custom Locus":
                locus = parsed_locus
            if title == "Variant Analysis":
                title = parsed_title
        else:
            raw_res = NC_ENGINE.analyze_fasta_entry(fasta_raw, wt_seq=None, flank_bp=flank_bp)
            s_raw = float(raw_res.get("score", 0.0))
            is_driver = raw_res.get("classification") in ("Driver Variant", "Loss-of-Function Variant")
            s_unified = s_raw if is_driver else 0.0
            return {
                **raw_res,
                "version": "NC-MutScan v3.3",
                "build_sha256": CORE_SHORT_SHA256,
                "build_verified": CORE_SHORT_SHA256 == EXPECTED_BUILD_PREFIX,
                "S_raw": s_raw,
                "S_unified": s_unified,
                "is_driver": is_driver,
            }

    if not wt_seq or not mut_seq:
        return None

    if raw_flank is None and len(wt_seq) == len(mut_seq) and len(wt_seq) >= 6:
        diffs = [i for i, (a, b) in enumerate(zip(wt_seq.upper(), mut_seq.upper())) if a != b]
        if diffs:
            flank_bp = diffs[0]

    raw_res = NC_ENGINE.analyze_variant_delta(
        wt_seq=wt_seq,
        mut_seq=mut_seq,
        flank_bp=flank_bp,
        locus=locus,
        title=title,
    )
    delta_pwm = float(raw_res.get("delta_pwm", 0.0))
    s_raw = round(abs(delta_pwm), 4)
    is_driver = raw_res.get("classification") in ("Driver Variant", "Loss-of-Function Variant")
    s_unified = s_raw if is_driver else 0.0

    return {
        **raw_res,
        "version": "NC-MutScan v3.3",
        "build_sha256": CORE_SHORT_SHA256,
        "build_verified": CORE_SHORT_SHA256 == EXPECTED_BUILD_PREFIX,
        "S_raw": s_raw,
        "S_unified": s_unified,
        "is_driver": is_driver,
    }


class NCMutScanRequestHandler(http.server.SimpleHTTPRequestHandler):
    def _set_headers(self, status=200, content_type="application/json"):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type, Authorization")
        self.end_headers()

    def do_OPTIONS(self):
        self._set_headers(200)

    def do_GET(self):
        path = urlparse(self.path).path
        if path in ("/", "/health", "/api/health"):
            self._set_headers(200)
            resp = {
                "status": "online",
                "engine": "NC-MutScan v3.3",
                "version": "NC-MutScan v3.3",
                "build_sha256": CORE_SHORT_SHA256,
                "full_sha256": CORE_FULL_SHA256,
                "build_verified": CORE_SHORT_SHA256 == EXPECTED_BUILD_PREFIX,
                "endpoints": ["/health", "/predict", "/api/analyze"],
            }
            self.wfile.write(json.dumps(resp, indent=2).encode("utf-8"))
            return

        if path == "/index.html":
            index_file = os.path.join(BASE_DIR, "index.html")
            if os.path.exists(index_file):
                self._set_headers(200, content_type="text/html; charset=utf-8")
                with open(index_file, "rb") as f:
                    self.wfile.write(f.read())
                return

        super().do_GET()

    def do_POST(self):
        path = urlparse(self.path).path
        try:
            content_length = int(self.headers.get("Content-Length", 0))
            if content_length > 2 * 1024 * 1024:
                self._set_headers(413)
                self.wfile.write(json.dumps({"error": "Payload exceeds 2 MB limit."}).encode("utf-8"))
                return

            post_data = self.rfile.read(content_length).decode("utf-8") if content_length > 0 else "{}"
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
                self._set_headers(200)
                self.wfile.write(json.dumps(res, indent=2).encode("utf-8"))
                return

            if path in ("/", "/predict", "/api/predict", "/api/analyze", "/analyze"):
                result = run_v33_prediction(payload)
                if result is None:
                    self._set_headers(400)
                    self.wfile.write(
                        json.dumps({
                            "error": "Provide 'wt_seq' and 'mut_seq' (or 'fasta' with header >chr:pos | Name (REF>ALT)).",
                            "version": "NC-MutScan v3.3",
                            "build_sha256": CORE_SHORT_SHA256,
                        }).encode("utf-8")
                    )
                    return

                self._set_headers(200)
                self.wfile.write(json.dumps(result, indent=2).encode("utf-8"))
                return

            self._set_headers(404)
            self.wfile.write(json.dumps({"error": f"Endpoint not found: {path}"}).encode("utf-8"))

        except json.JSONDecodeError:
            self._set_headers(400)
            self.wfile.write(json.dumps({"error": "Invalid JSON payload."}).encode("utf-8"))
        except Exception as e:
            self._set_headers(500)
            self.wfile.write(json.dumps({"error": str(e)}).encode("utf-8"))



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


def run():
    server_address = ("0.0.0.0", PORT)
    httpd = http.server.ThreadingHTTPServer(server_address, NCMutScanRequestHandler)
    print("==========================================================================", flush=True)
    print(f" NC-MutScan v3.3 API Server | build sha256:{CORE_SHORT_SHA256}", flush=True)
    print(f" Listening on http://0.0.0.0:{PORT} (/health, /predict, /api/analyze)", flush=True)
    print("==========================================================================", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        httpd.server_close()


if __name__ == "__main__":
    run()
