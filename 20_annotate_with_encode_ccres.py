import gzip
import json
import os
from collections import Counter, defaultdict

MODEL_PATH = "./trained_server_model/noncoding_driver_db.json"
UPDATED_MODEL_PATH = "./trained_server_model/noncoding_driver_db.json"
RESULTS_DIR = "./results"

# Epigenomic functional weight multipliers for ENCODE cCRE signatures
CCRE_WEIGHTS = {
    "PLS": 2.5,        # Promoter-like signature
    "pELS": 2.0,       # Proximal Enhancer-like signature
    "dELS": 2.0,       # Distal Enhancer-like signature
    "CTCF-only": 1.5,  # Insulator / Chromatin boundary
    "CTCF": 1.5,
    "DNase-only": 1.2, # Open chromatin region
    "DNase": 1.2
}

print("[*] Locating ENCODE cCRE dataset...", flush=True)

target_bed_path = None
search_locations = [
    RESULTS_DIR, 
    ".", 
    "./ENCFF924IMH.bed", 
    os.path.join(RESULTS_DIR, "ENCFF924IMH.bed")
]

for loc in search_locations:
    if os.path.isfile(loc) and (loc.endswith('.bed') or loc.endswith('.bed.gz')):
        target_bed_path = loc
        break
    elif os.path.isdir(loc):
        for root, dirs, files in os.walk(loc):
            for f in files:
                if f.endswith('.bed') or f.endswith('.bed.gz'):
                    target_bed_path = os.path.join(root, f)
                    break
            if target_bed_path:
                break
    if target_bed_path:
        break

if not target_bed_path:
    print("[!] ERROR: Could not locate .bed or .bed.gz file!", flush=True)
    exit(1)

print(f"[+] Found cCRE dataset: {target_bed_path}", flush=True)

# ---------------------------------------------------------
# 1. LOAD TRAINED KNOWLEDGE BASE
# ---------------------------------------------------------
if not os.path.exists(MODEL_PATH):
    print(f"[!] ERROR: Knowledge base not found at {MODEL_PATH}.", flush=True)
    exit(1)

with open(MODEL_PATH, "r", encoding="utf-8") as f:
    knowledge_base = json.load(f)

# ---------------------------------------------------------
# 2. 10KB GENOMIC INTERVAL BUCKETING (FAST O(1) LOOKUP)
# ---------------------------------------------------------
BUCKET_SIZE = 10000
chrom_buckets = defaultdict(list)

for region_id, item in knowledge_base.items():
    c = item["chrom"]
    s = item["start"]
    e = item["end"]
    start_b = s // BUCKET_SIZE
    end_b = e // BUCKET_SIZE
    for b in range(start_b, end_b + 1):
        chrom_buckets[(c, b)].append((s, e, region_id))

print(f"[+] Indexed {len(knowledge_base):,} non-coding loci into 10kb interval buckets.", flush=True)

# ---------------------------------------------------------
# 3. FAST BUCKETED SCAN OVER 1M+ ENCODE cCREs
# ---------------------------------------------------------
print("[*] Scanning 1M+ ENCODE cCREs using Fast Bucket Index...", flush=True)

ccre_matches = defaultdict(list)
ccre_type_counts = Counter()
total_ccres_scanned = 0

open_fn = gzip.open if target_bed_path.endswith('.gz') else open

with open_fn(target_bed_path, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("#") or line.startswith("track") or not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue

        chrom = parts[0]
        try:
            c_start = int(parts[1])
            c_end = int(parts[2])
        except ValueError:
            continue

        full_line_str = " ".join(parts[3:]) if len(parts) > 3 else ""
        ccre_type = "Low-DNase"
        
        for label_key in ["PLS", "pELS", "dELS", "CTCF-only", "CTCF", "DNase-only", "DNase"]:
            if label_key in full_line_str:
                ccre_type = label_key.replace("-only", "")
                break

        total_ccres_scanned += 1

        # Only query loci in 10kb buckets touched by this cCRE
        start_b = c_start // BUCKET_SIZE
        end_b = c_end // BUCKET_SIZE
        
        seen_in_line = set()
        for b in range(start_b, end_b + 1):
            for l_start, l_end, region_id in chrom_buckets.get((chrom, b), []):
                if region_id not in seen_in_line:
                    seen_in_line.add(region_id)
                    if max(c_start, l_start) < min(c_end, l_end):
                        ccre_matches[region_id].append(ccre_type)
                        ccre_type_counts[ccre_type] += 1

print(f"[+] Scanned {total_ccres_scanned:,} ENCODE cCRE regions.", flush=True)

# ---------------------------------------------------------
# 4. RE-SCORE HOTSPOTS WITH ENCODE WEIGHTS
# ---------------------------------------------------------
print("[*] Re-calculating pathogenicity scores...", flush=True)

updated_driver_count = 0

for region_id, item in knowledge_base.items():
    matched_ccres = list(set(ccre_matches.get(region_id, [])))
    item["encode_ccres"] = matched_ccres

    ccre_bonus = 1.0
    if matched_ccres:
        ccre_bonus = max([CCRE_WEIGHTS.get(ct, 1.0) for ct in matched_ccres])

    base_score = item["score"]
    revised_score = round(base_score * ccre_bonus, 2)
    item["score"] = revised_score

    is_driver = revised_score >= 4.0 or (item["recurrence"] >= 2 and len(matched_ccres) > 0 and ccre_bonus >= 2.0)
    item["is_driver"] = is_driver

    if is_driver:
        updated_driver_count += 1

with open(UPDATED_MODEL_PATH, "w", encoding="utf-8") as f:
    json.dump(knowledge_base, f, indent=2)

print("\n==========================================================================")
print(" ENCODE cCRE INTEGRATION COMPLETE")
print("==========================================================================")
print(f"[+] Total Hotspots Annotated with ENCODE cCREs : {len(ccre_matches):,}")
print(f"[+] High-Confidence Driver Loci Count          : {updated_driver_count:,}")
print(f"[+] Updated Knowledge Base Saved To            : {UPDATED_MODEL_PATH}")
print("\n--- Overlapping ENCODE Regulatory Signature Hits ---")
for ctype, count in ccre_type_counts.most_common():
    print(f"  • {ctype:<12}: {count:>5,} locus matches")
print("==========================================================================")