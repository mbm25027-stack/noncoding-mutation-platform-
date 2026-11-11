import csv
import json
import os

INPUT_HOTSPOTS = "./results/noncoding_hotspots_summary.csv"
MODEL_DIR = "./trained_server_model"
OUTPUT_MODEL_JSON = os.path.join(MODEL_DIR, "noncoding_driver_db.json")

os.makedirs(MODEL_DIR, exist_ok=True)

print("[*] Annotating non-coding hotspots and training server rules...", flush=True)

# Weight factors for variant types based on biological impact severity
IMPACT_WEIGHTS = {
    "Splice_Site": 2.5,
    "5'UTR": 2.2,
    "5'Flank": 2.0,  # Promoter regions
    "Splice_Region": 1.8,
    "3'UTR": 1.5,
    "RNA": 1.3,
    "3'Flank": 1.2,
    "Intron": 0.8
}

knowledge_base = {}

total_processed = 0
high_confidence_drivers = 0

with open(INPUT_HOTSPOTS, 'r', encoding='utf-8', errors='ignore') as f:
    reader = csv.DictReader(f)
    for row in reader:
        total_processed += 1
        chrom = row.get('Chromosome', '')
        try:
            win_start = int(row.get('Window_Start', 0))
            win_end = int(row.get('Window_End', 0))
            recurrence = int(row.get('Recurrence_Count', 0))
            samples = int(row.get('Sample_Count', 0))
        except (ValueError, TypeError):
            continue

        gene = row.get('Gene', 'Unknown')
        v_types = row.get('Variant_Types', '').split(';')
        cancer_types = row.get('Cancer_Types', '').split(';')

        # Calculate functional score
        max_impact = max([IMPACT_WEIGHTS.get(vt, 1.0) for vt in v_types]) if v_types else 1.0
        score = round((recurrence * 0.6 + samples * 0.4) * max_impact, 2)

        is_driver = score >= 4.0 or (recurrence >= 3 and max_impact >= 1.8)

        if is_driver:
            high_confidence_drivers += 1

        region_id = f"{chrom}:{win_start}-{win_end}"
        knowledge_base[region_id] = {
            "chrom": chrom,
            "start": win_start,
            "end": win_end,
            "gene": gene,
            "recurrence": recurrence,
            "samples": samples,
            "score": score,
            "is_driver": is_driver,
            "variant_types": v_types,
            "cancer_types": cancer_types
        }

# Save trained database for server lookup
with open(OUTPUT_MODEL_JSON, 'w', encoding='utf-8') as f:
    json.dump(knowledge_base, f, indent=2)

print(f"[+] SUCCESS: Processed {total_processed} regions.", flush=True)
print(f"[+] Identified {high_confidence_drivers} high-confidence non-coding driver loci.", flush=True)
print(f"[+] Non-coding Knowledge Base compiled at: {OUTPUT_MODEL_JSON}", flush=True)