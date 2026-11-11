import http.server
import json
import os
import sys
from urllib.parse import urlparse

# Import sequence engine with graceful fallback
align_and_analyze_sequence = None
try:
    from sequence_aligner_engine import align_and_analyze_sequence
except ImportError:
    try:
        import importlib
        seq_mod = importlib.import_module("sequence_aligner_engine")
        align_and_analyze_sequence = seq_mod.align_and_analyze_sequence
    except Exception as e:
        print(f"[!] Warning: sequence_aligner_engine could not be loaded: {e}", flush=True)

# Import Epigenomics engine
try:
    import epigenomics_engine
except ImportError:
    epigenomics_engine = None

# Import Structural Variant engine
try:
    import sv_engine
except ImportError:
    sv_engine = None

PORT = int(os.environ.get("PORT", 8000))
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(BASE_DIR, "trained_server_model", "noncoding_driver_db.json")

NONCODING_DB = {}

if os.path.exists(MODEL_PATH):
    try:
        with open(MODEL_PATH, "r", encoding="utf-8") as f:
            NONCODING_DB = json.load(f)
        print(f"[+] Loaded {len(NONCODING_DB):,} non-coding loci into memory.", flush=True)
    except Exception as e:
        print(f"[!] ERROR loading model database: {e}", flush=True)


def normalize_chrom(chrom):
    if not chrom:
        return ""
    c_str = str(chrom).strip()
    return c_str if c_str.lower().startswith("chr") else f"chr{c_str}"


def analyze_mutation(chrom, pos, gene=None):
    norm_chrom = normalize_chrom(chrom)
    try:
        pos_int = int(pos)
    except (ValueError, TypeError):
        return {"error": "Invalid genomic position integer."}

    # Fetch epigenomic context (open chromatin / cCRE / CTCF)
    epi_context = None
    if epigenomics_engine and hasattr(epigenomics_engine, "query_locus_context"):
        try:
            epi_context = epigenomics_engine.query_locus_context(norm_chrom, pos_int)
        except Exception as e:
            print(f"[!] Epigenomic query error: {e}", flush=True)

    matched_hit = None
    for region_id, data in NONCODING_DB.items():
        if data.get("chrom") == norm_chrom:
            if data.get("start") <= pos_int <= data.get("end"):
                matched_hit = data
                break

    if matched_hit:
        is_driver = matched_hit.get("is_driver", False)
        base_score = matched_hit.get("score", 0.0)
        multiplier = epi_context.get("multiplier", 1.0) if isinstance(epi_context, dict) else 1.0
        final_score = round(base_score * multiplier, 4)

        return {
            "status": "Pathogenic / Non-Coding Driver" if is_driver else "Low-Impact Non-Coding Locus",
            "is_driver": is_driver,
            "score": final_score,
            "base_score": base_score,
            "epigenomic_multiplier": multiplier,
            "epigenomic_context": epi_context,
            "gene": matched_hit.get("gene", "Unknown"),
            "recurrence": matched_hit.get("recurrence", 0),
            "samples": matched_hit.get("samples", 0),
            "cancer_types": matched_hit.get("cancer_types", []),
            "variant_types": matched_hit.get("variant_types", []),
            "encode_ccres": matched_hit.get("encode_ccres", []),
            "region": f"{matched_hit.get('chrom')}:{matched_hit.get('start')}-{matched_hit.get('end')}",
            "message": "Locus matched high-confidence non-coding driver region."
        }

    return {
        "status": "Benign / Safe Locus",
        "is_driver": False,
        "score": 0.0,
        "epigenomic_context": epi_context,
        "gene": gene or "Unknown",
        "recurrence": 0,
        "samples": 0,
        "cancer_types": [],
        "variant_types": [],
        "encode_ccres": [],
        "region": f"{norm_chrom}:{pos_int}",
        "message": "No recurrent non-coding driver hotspot matched at this position."
    }


class NCDriverRequestHandler(http.server.SimpleHTTPRequestHandler):

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
        parsed_url = urlparse(self.path)
        path = parsed_url.path

        if path in ["/health", "/api/health"]:
            self._set_headers(200)
            is_epi_loaded = getattr(epigenomics_engine, "IS_LOADED", False) if epigenomics_engine else False
            response = {
                "status": "online",
                "engine": "NC-DRIVER PAN-CANCER DIAGNOSTIC SERVER",
                "version": "3.1-EPIGENOMIC-SV-MULTITHREADED",
                "loaded_driver_regions": len(NONCODING_DB),
                "epigenomics_indexed": is_epi_loaded
            }
            self.wfile.write(json.dumps(response, indent=2).encode("utf-8"))
            return

        if path in ["/", "/index.html"]:
            index_file = os.path.join(BASE_DIR, "index.html")
            if os.path.exists(index_file):
                self._set_headers(200, content_type="text/html")
                with open(index_file, "rb") as f:
                    self.wfile.write(f.read())
                return

        super().do_GET()

    def do_POST(self):
        parsed_url = urlparse(self.path)
        path = parsed_url.path

        # 1. Structural Variant API Route
        if path in ["/api/analyze_sv", "/analyze_sv"]:
            try:
                content_length = int(self.headers.get("Content-Length", 0))
                if content_length > 2 * 1024 * 1024:
                    self._set_headers(413)
                    self.wfile.write(json.dumps({"error": "Payload exceeds 2 MB limit."}).encode("utf-8"))
                    return

                post_data = self.rfile.read(content_length).decode("utf-8")
                payload = json.loads(post_data) if post_data else {}

                sv_type = payload.get("sv_type", "TRANSLOCATION")
                chrom_a = payload.get("chrom_a", "")
                pos_a = payload.get("pos_a", 0)
                chrom_b = payload.get("chrom_b", "")
                pos_b = payload.get("pos_b", 0)

                if not chrom_a or not pos_a:
                    self._set_headers(400)
                    self.wfile.write(json.dumps({"error": "Missing Breakpoint A coordinates (chrom_a, pos_a)."}).encode("utf-8"))
                    return

                if sv_engine and hasattr(sv_engine, "analyze_structural_variant"):
                    result = sv_engine.analyze_structural_variant(
                        sv_type, chrom_a, int(pos_a), chrom_b, int(pos_b)
                    )
                else:
                    result = {"error": "Structural Variant engine (sv_engine.py) not available."}

                self._set_headers(200)
                self.wfile.write(json.dumps(result, indent=2).encode("utf-8"))
                return
            except json.JSONDecodeError:
                self._set_headers(400)
                self.wfile.write(json.dumps({"error": "Invalid JSON payload."}).encode("utf-8"))
            except Exception as e:
                self._set_headers(500)
                self.wfile.write(json.dumps({"error": str(e)}).encode("utf-8"))
            return

        # 2. Standard SNV / Sequence API Route
        if path in ["/", "/api/analyze", "/analyze"]:
            try:
                content_length = int(self.headers.get("Content-Length", 0))

                if content_length > 2 * 1024 * 1024:
                    self._set_headers(413)
                    self.wfile.write(json.dumps({"error": "Payload exceeds 2 MB limit."}).encode("utf-8"))
                    return

                post_data = self.rfile.read(content_length).decode("utf-8")
                payload = json.loads(post_data) if post_data else {}

                raw_sequence = payload.get("sequence", "").strip()

                # Dual Mode A: Sequence analysis
                if raw_sequence:
                    if callable(align_and_analyze_sequence):
                        result = align_and_analyze_sequence(raw_sequence, NONCODING_DB)
                    else:
                        result = {"error": "Sequence alignment engine unavailable on server."}

                    self._set_headers(200)
                    self.wfile.write(json.dumps(result, indent=2).encode("utf-8"))
                    return

                # Dual Mode B: Coordinate analysis
                chrom = payload.get("chrom", "")
                pos = payload.get("pos", 0)
                gene = payload.get("gene", "")

                if not chrom or not pos:
                    self._set_headers(400)
                    self.wfile.write(json.dumps({"error": "Provide either 'sequence' OR ('chrom' and 'pos')."}).encode("utf-8"))
                    return

                result = analyze_mutation(chrom, pos, gene)
                self._set_headers(200)
                self.wfile.write(json.dumps(result, indent=2).encode("utf-8"))

            except json.JSONDecodeError:
                self._set_headers(400)
                self.wfile.write(json.dumps({"error": "Invalid JSON payload."}).encode("utf-8"))
            except Exception as e:
                self._set_headers(500)
                self.wfile.write(json.dumps({"error": str(e)}).encode("utf-8"))
        else:
            self._set_headers(404)
            self.wfile.write(json.dumps({"error": "Endpoint not found."}).encode("utf-8"))


def run():
    server_address = ("0.0.0.0", PORT)
    httpd = http.server.ThreadingHTTPServer(server_address, NCDriverRequestHandler)
    
    is_epi_active = getattr(epigenomics_engine, "IS_LOADED", False) if epigenomics_engine else False

    print("==========================================================================", flush=True)
    print(" NC-DRIVER ENGINE | MULTITHREADED PAN-CANCER SERVER v3.1", flush=True)
    print(f" [*] Server running on: http://0.0.0.0:{PORT}", flush=True)
    print(f" [*] Loaded driver regions: {len(NONCODING_DB):,}", flush=True)
    print(f" [*] ENCODE Epigenomics Engine: {'ACTIVE' if is_epi_active else 'OFFLINE'}", flush=True)
    print(f" [*] Structural Variant Engine: {'ACTIVE' if sv_engine else 'OFFLINE'}", flush=True)
    print("==========================================================================", flush=True)
    
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n[-] Shutting down server gracefully.", flush=True)
        httpd.server_close()


if __name__ == "__main__":
    run()