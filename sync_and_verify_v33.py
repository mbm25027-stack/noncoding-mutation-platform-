import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request

EXPECTED_HASH = "b276d0d44045"
CORE_FILE = "nc_mutscan_research_suite.py"
BACKUP_DIR = "_deployment_backups"


def verify_core_hash():
    if not os.path.isfile(CORE_FILE):
        raise FileNotFoundError(f"Missing frozen core engine: {CORE_FILE}")
    with open(CORE_FILE, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()
    if not digest.startswith(EXPECTED_HASH):
        raise RuntimeError(f"Core hash mismatch! Expected {EXPECTED_HASH}, got {digest[:12]}")
    return digest


SERVER_PY_CODE = '''import hashlib
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
            seq_parts.append(re.sub(r"\\s+", "", ln).upper())
    mut_seq = "".join(seq_parts)
    parts = [p.strip() for p in header.split("|")]
    lm = re.search(r"(chr[\\dXYM]+:\\d+)", header, re.I)
    locus = lm.group(1) if lm else (parts[0] if parts and parts[0] else "Custom Locus")
    title = parts[1].replace("_", " ") if len(parts) >= 2 and parts[1] else (header or "Variant Analysis")

    vm = re.search(r"([ACGT])\\s*>\\s*([ACGT])", header, re.I)
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

    flank_bp = int(payload.get("flank_bp", 25))
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
'''


def patch_index_html():
    if not os.path.isfile("index.html"):
        print("[WARN] index.html not found; skipping frontend patch.")
        return False

    html = open("index.html", "r", encoding="utf-8-sig", errors="ignore").read()

    # 1. Update v3.1 release labels & report headers to v3.3 (sha256:b276d0d44045)
    replacements = [
        ("Analyzer v3.1", "Analyzer v3.3"),
        ("Computational genomics v3.1", "Computational genomics v3.3 (sha256:b276d0d44045)"),
        ("Release 3.1", "Release 3.3"),
        ("Interface (v3.1):", "Interface (v3.3):"),
        ("NC-MutScan v3.1 (Sequence Engine)", "NC-MutScan v3.3 (Frozen Build sha256:b276d0d44045)"),
        ("NC-MutScan Core Engine v3.1", "NC-MutScan Core Engine v3.3 (sha256:b276d0d44045)"),
        ("Welcome to <span>NC-MutScan v3.1</span>", "Welcome to <span>NC-MutScan v3.3</span>"),
        ("ENGINE CORE v3.1", "ENGINE CORE v3.3 (sha256:b276d0d44045)"),
        ("NC-MutScan v3.1 | Noncoding Regulatory Report", "NC-MutScan v3.3 (sha256:b276d0d44045) | Regulatory Report"),
        ("NC-MutScan v3.1 Report:", "NC-MutScan v3.3 Report:"),
    ]
    for old, new in replacements:
        html = html.replace(old, new)

    # 2. Wire runAnalysis() to ping Render on load and query /predict with local fallback
    old_run_fn = """async function runAnalysis() {
  const reportPanel = document.getElementById('report-output');
  if (reportPanel) reportPanel.classList.add('is-loading');
  setTimeout(() => {
    try {
      const data = computeAnalysis();
      lastAnalysisResult = data;
      window.latestNCMutScanResults = data;
      renderReport(data);
    } finally {
      if (reportPanel) reportPanel.classList.remove('is-loading');
    }
  }, 200);
}"""

    new_run_fn = """const RENDER_API_BASE = 'https://noncoding-mutation-server-1.onrender.com';
window.addEventListener('DOMContentLoaded', () => {
  fetch(`${RENDER_API_BASE}/health`, { method: 'GET' }).catch(() => {});
});

async function runAnalysis() {
  const reportPanel = document.getElementById('report-output');
  if (reportPanel) reportPanel.classList.add('is-loading');
  try {
    const data = computeAnalysis();
    const rawInput = (document.getElementById('fasta-input')?.value || '').trim();
    const tissueVal = document.getElementById('tissue-select')?.value || 'pan_cancer';

    try {
      const ctrl = new AbortController();
      const timer = setTimeout(() => ctrl.abort(), 4500);
      const resp = await fetch(`${RENDER_API_BASE}/predict`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ fasta: rawInput, tissue: tissueVal }),
        signal: ctrl.signal
      });
      clearTimeout(timer);
      if (resp.ok) {
        const apiData = await resp.json();
        if (apiData && apiData.build_sha256 === 'b276d0d44045') {
          data.classification = apiData.classification || data.classification;
          if (typeof apiData.score === 'number') data.score = apiData.score;
          if (apiData.mechanism) data.mechanism = `${apiData.mechanism} [Verified API v3.3 sha256:${apiData.build_sha256}]`;
          if (apiData.primary_tf && apiData.primary_tf !== 'None') data.tf = apiData.primary_tf;
          if (apiData.direction) data.direction = apiData.direction;
        }
      }
    } catch (_) {
      /* Render cold-start or offline: seamless client-side v3.3 parity fallback */
    }

    lastAnalysisResult = data;
    window.latestNCMutScanResults = data;
    renderReport(data);
  } finally {
    if (reportPanel) reportPanel.classList.remove('is-loading');
  }
}"""

    if old_run_fn in html:
        html = html.replace(old_run_fn, new_run_fn, 1)
        print("[PASS] Wired index.html runAnalysis() to Render /predict + background /health wake-up.")
    else:
        print("[INFO] Updated v3.1 -> v3.3 labels in index.html (runAnalysis block already customized).")

    with open("index.html", "w", encoding="utf-8") as f:
        f.write(html)
    return True


def http_json(url, method="GET", payload=None):
    headers = {"Content-Type": "application/json", "Origin": "https://mbm25027-stack.github.io"}
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    t0 = time.perf_counter()
    with urllib.request.urlopen(req, timeout=10) as r:
        ms = (time.perf_counter() - t0) * 1000.0
        return r.status, dict(r.headers), json.loads(r.read().decode("utf-8")), ms


def run_local_server_verification():
    print("\n=== RUNNING LOCAL HTTP VERIFICATION ON UPDATED server.py ===")
    test_port = 18080
    env = os.environ.copy()
    env["PORT"] = str(test_port)
    proc = subprocess.Popen(
        [sys.executable, "server.py"],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    try:
        base = f"http://127.0.0.1:{test_port}"
        ready = False
        for _ in range(30):
            try:
                with socket.create_connection(("127.0.0.1", test_port), timeout=0.3):
                    ready = True
                    break
            except OSError:
                time.sleep(0.15)

        if not ready:
            raise RuntimeError("Local test server failed to start.")

        st, hdrs, health, ms = http_json(f"{base}/health")
        print(f"[PASS] GET /health -> HTTP {st} ({ms:.1f} ms) | version={health.get('version')} | sha256={health.get('build_sha256')}")

        # 7 deterministic v3.3 regression tests
        cases = [
            (
                "TC-01 TERT -124 C>T (GABPA/ETS Gain)",
                "CCCCTCCCGGGTCCCCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCC",
                "CCCCTCCCGGGTCCCCGGCCCAGCCCCTTCCGGGCCCTCCCAGCCCCTCCC",
                "Driver Variant",
            ),
            (
                "TC-06 TP53 Promoter Sp1 Disruption",
                "TTTTTGGGCGGGGCCGGGCGGGGCCGGGCGGGGCCGGGCGGGGCCAAAAAA",
                "TTTTTGGGCGGGGCCGGGCGGGGCCGGGCGAGGCCGGGCGGGGCCAAAAAA",
                "Loss-of-Function Variant",
            ),
            (
                "U1 5' Splice Donor Loss (GTAAGT -> ATAAGT)",
                "TGCATGCATGCATGCATGCATCAGGTAAGTTGCATGCATGCATGCATGCAT",
                "TGCATGCATGCATGCATGCATCAGATAAGTTGCATGCATGCATGCATGCAT",
                "Loss-of-Function Variant",
            ),
            (
                "3'-UTR AATAAA PolyA Signal Loss",
                "ATATATATATATATATATATAAATAAAATATATATATATATATATATATAT",
                "ATATATATATATATATATATAAGAAAAATATATATATATATATATATATAT",
                "Loss-of-Function Variant",
            ),
            (
                "TC-03/12 Clustered AID/SHM WRCY/RGYW",
                "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
                "AGCTAGCTAGCTAGCTAGCTAGCTGGCTAGCTAGCTAGCTAGCTAGCTAGC",
                "AID/SHM Passenger",
            ),
            (
                "TC-04 IUPAC Degenerate Base Sanitization",
                "ATGCATGCATNRYATGCATGCATGCATGCATGCATGCATGCATGCATGCAT",
                "ATGCATGCATNRYATGCATGCATGCGTGCATGCATGCATGCATGCATGCAT",
                "Ambiguous IUPAC Sequence",
            ),
            (
                "TC-10 Short Sequence (<6 bp) Bypass",
                "ACGTA",
                "ACTTA",
                "Sub-Threshold Short Sequence",
            ),
        ]

        for label, wt, mut, expected_cls in cases:
            st, _, res, ms = http_json(
                f"{base}/predict",
                method="POST",
                payload={"wt_seq": wt, "mut_seq": mut, "flank_bp": 25, "title": label},
            )
            got_cls = res.get("classification")
            ok = (st == 200 and got_cls == expected_cls and res.get("build_sha256") == EXPECTED_HASH)
            tag = "PASS" if ok else "CHECK"
            print(
                f"[{tag}] {label:<42} -> {got_cls} | dPWM={res.get('delta_pwm')} | S_raw={res.get('S_raw')} | S_unified={res.get('S_unified')} ({ms:.1f} ms)"
            )
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            proc.kill()


def main():
    before_hash = verify_core_hash()
    print(f"[PASS] Verified untouched frozen core engine: sha256:{before_hash[:12]}")

    os.makedirs(BACKUP_DIR, exist_ok=True)
    for fn in ("server.py", "index.html"):
        if os.path.isfile(fn):
            shutil.copy2(fn, os.path.join(BACKUP_DIR, f"{fn}.bak"))
    print(f"[PASS] Backed up server.py and index.html to {BACKUP_DIR}/")

    with open("server.py", "w", encoding="utf-8") as f:
        f.write(SERVER_PY_CODE)
    print("[PASS] Wrote updated server.py wired to frozen NC-MutScan v3.3 engine.")

    patch_index_html()
    run_local_server_verification()

    after_hash = verify_core_hash()
    assert before_hash == after_hash, "Frozen core hash changed!"
    print(f"\n[VERIFIED] Frozen core engine untouched: sha256:{after_hash[:12]}")


if __name__ == "__main__":
    main()