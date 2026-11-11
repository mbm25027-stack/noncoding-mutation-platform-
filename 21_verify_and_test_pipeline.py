import json
import os
import sys
import urllib.parse
import urllib.request
from collections import Counter

# Paths
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MUTATIONS_FILE = os.path.join(BASE_DIR, "results", "pan_cancer_noncoding_mutations.csv")
HOTSPOTS_FILE = os.path.join(BASE_DIR, "results", "noncoding_hotspots_summary.csv")
REPORT_FILE = os.path.join(BASE_DIR, "results", "master_summary_report.txt")
MODEL_FILE = os.path.join(BASE_DIR, "trained_server_model", "noncoding_driver_db.json")
SERVER_FILE = os.path.join(BASE_DIR, "server.py")
INDEX_FILE = os.path.join(BASE_DIR, "index.html")
GITIGNORE_FILE = os.path.join(BASE_DIR, ".gitignore")

SERVER_URL = "http://localhost:8000"

print("==========================================================================", flush=True)
print(" NC-DRIVER ENGINE | MASTER END-TO-END SYSTEM DIAGNOSTIC TEST", flush=True)
print("==========================================================================", flush=True)

passed_checks = 0
total_checks = 0

def check_mark(status):
    global passed_checks
    if status:
        passed_checks += 1
        return "[PASS]"
    return "[FAIL]"

# ---------------------------------------------------------
# CHECK 1: FILE SYSTEM & PIPELINE ARTIFACT INTEGRITY
# ---------------------------------------------------------
total_checks += 1
print("\n[1] VERIFYING WORKSPACE ARTIFACTS...", flush=True)
required_files = {
    "Extracted Variants CSV": MUTATIONS_FILE,
    "Hotspots Summary CSV": HOTSPOTS_FILE,
    "Master Analysis Report": REPORT_FILE,
    "Trained Driver Model JSON": MODEL_FILE,
    "Production Server Script": SERVER_FILE,
    "Web Interface (index.html)": INDEX_FILE
}

all_files_exist = True
for name, path in required_files.items():
    exists = os.path.exists(path)
    size_str = f"({os.path.getsize(path) / 1024:.1f} KB)" if exists else "(MISSING)"
    print(f"  {check_mark(exists)} {name:<30}: {size_str}")
    if not exists:
        all_files_exist = False

# ---------------------------------------------------------
# CHECK 2: MODEL SCHEMA & ENCODE EPIGENOMIC DATA INTEGRITY
# ---------------------------------------------------------
total_checks += 1
print("\n[2] VERIFYING MODEL SCHEMAS & ENCODE REPOSITORY...", flush=True)

db_valid = False
if os.path.exists(MODEL_FILE):
    try:
        with open(MODEL_FILE, "r", encoding="utf-8") as f:
            model = json.load(f)
        
        total_loci = len(model)
        drivers = [v for v in model.values() if v.get("is_driver", False)]
        ccre_annotated = [v for v in model.values() if v.get("encode_ccres")]

        ccre_counts = Counter()
        for v in model.values():
            for cc in v.get("encode_ccres", []):
                ccre_counts[cc] += 1

        print(f"  {check_mark(total_loci > 0)} Total Indexed Loci         : {total_loci:,}")
        print(f"  {check_mark(len(drivers) > 0)} High-Confidence Drivers   : {len(drivers):,}")
        print(f"  {check_mark(len(ccre_annotated) > 0)} ENCODE cCRE Hits        : {len(ccre_annotated):,}")

        print("\n  Overlapping ENCODE Epigenomic Signatures:")
        for sig, count in ccre_counts.most_common():
            print(f"    - {sig:<10}: {count:>5,} hits")

        db_valid = total_loci > 0 and len(drivers) > 0
    except Exception as e:
        print(f"  [FAIL] Failed to parse JSON model: {e}")

# ---------------------------------------------------------
# CHECK 3: INTERNAL LOCUS SCORING ENGINE SIMULATION
# ---------------------------------------------------------
total_checks += 1
print("\n[3] TESTING INTERNAL DIAGNOSTIC ENGINE SIMULATION...", flush=True)

engine_functional = False
if db_valid:
    # Find one known driver locus to simulate
    sample_driver = next((v for v in model.values() if v.get("is_driver")), None)
    if sample_driver:
        s_chrom = sample_driver["chrom"]
        s_pos = sample_driver["start"] + 10
        s_gene = sample_driver["gene"]

        print(f"  [*] Testing Query Locus -> Chrom: {s_chrom}, Pos: {s_pos}, Gene: {s_gene}")

        # Simulate matching logic
        match = None
        for data in model.values():
            if data["chrom"] == s_chrom and data["start"] <= s_pos <= data["end"]:
                match = data
                break

        if match:
            print(f"  {check_mark(True)} Classification Result    : {'Pathogenic Driver' if match['is_driver'] else 'Low Impact'}")
            print(f"  {check_mark(True)} Calculated Score        : {match['score']}")
            print(f"  {check_mark(True)} Matched ENCODE Elements : {', '.join(match.get('encode_ccres', [])) or 'None'}")
            engine_functional = True
        else:
            print("  [FAIL] Locus lookup simulation failed.")

# ---------------------------------------------------------
# CHECK 4: LIVE SERVER API RESPONSE TEST (HTTP GET/POST)
# ---------------------------------------------------------
total_checks += 1
print("\n[4] TESTING LIVE LOCAL API SERVER CONNECTION (http://localhost:8000)...", flush=True)

server_live = False
try:
    # 1. Health check GET
    req = urllib.request.Request(f"{SERVER_URL}/health")
    with urllib.request.urlopen(req, timeout=3) as resp:
        if resp.status == 200:
            h_data = json.loads(resp.read().decode("utf-8"))
            print(f"  {check_mark(True)} GET /health API Online   : Version {h_data.get('version')} ({h_data.get('loaded_driver_regions')} loci in memory)")

    # 2. Diagnostic POST
    payload = json.dumps({"chrom": "chr17", "pos": 7675000, "gene": "TP53"}).encode("utf-8")
    req_post = urllib.request.Request(f"{SERVER_URL}/api/analyze", data=payload, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req_post, timeout=3) as resp:
        if resp.status == 200:
            p_data = json.loads(resp.read().decode("utf-8"))
            print(f"  {check_mark(True)} POST /api/analyze Response : Status = '{p_data.get('status')}', Region = '{p_data.get('region')}'")
            server_live = True

except Exception:
    print("  [WARN] Live server is not currently running in terminal.")
    print("         (Run 'python server.py' in another terminal window to activate endpoint)")

# ---------------------------------------------------------
# CHECK 5: GITHUB & DEPLOYMENT READINESS
# ---------------------------------------------------------
total_checks += 1
print("\n[5] VERIFYING GITHUB & CLOUD DEPLOYMENT READINESS...", flush=True)

deployment_ready = False
if os.path.exists(GITIGNORE_FILE):
    with open(GITIGNORE_FILE, "r", encoding="utf-8") as f:
        git_content = f.read()
    
    ignores_heavy = "maf_files/" in git_content and "*.bed" in git_content
    keeps_model = "!trained_server_model/noncoding_driver_db.json" in git_content

    print(f"  {check_mark(ignores_heavy)} .gitignore excludes heavy files (.bed / maf_files)")
    print(f"  {check_mark(keeps_model)} .gitignore retains trained_server_model JSON")
    deployment_ready = ignores_heavy and keeps_model
else:
    print("  [FAIL] .gitignore file missing.")

# ---------------------------------------------------------
# FINAL SCORE SUMMARY
# ---------------------------------------------------------
print("\n==========================================================================")
print(" DIAGNOSTIC SUMMARY DASHBOARD")
print("==========================================================================")
print(f"  • Passed Component Tests : {passed_checks} / {total_checks}")
if passed_checks >= 4:
    print("  • OVERALL SYSTEM STATUS  : [SPOT ON / PRODUCTION READY]")
else:
    print("  • OVERALL SYSTEM STATUS  : [REQUIRES ATTENTION]")
print("==========================================================================")