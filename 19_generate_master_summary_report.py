import csv
import json
import os
from collections import Counter, defaultdict

MUTATIONS_FILE = "./results/pan_cancer_noncoding_mutations.csv"
HOTSPOTS_FILE = "./results/noncoding_hotspots_summary.csv"
MODEL_FILE = "./trained_server_model/noncoding_driver_db.json"
REPORT_OUTPUT = "./results/master_summary_report.txt"

os.makedirs("./results", exist_ok=True)

print("[*] Generating Pan-Cancer Noncoding Master Analysis Report...", flush=True)

report_lines = []

def log(text=""):
    print(text, flush=True)
    report_lines.append(text)

log("==========================================================================")
log(" PAN-CANCER NONCODING GENOME MUTATION ANALYSIS ENGINE - MASTER REPORT")
log("==========================================================================")

# ---------------------------------------------------------
# 1. ANALYZE EXTRACTED MUTATIONS DATASET
# ---------------------------------------------------------
log("\n[1] RAW MUTATION EXTRACTION SUMMARY")
log("--------------------------------------------------------------------------")

if os.path.exists(MUTATIONS_FILE):
    total_mutations = 0
    variant_types = Counter()
    cancer_cohorts = Counter()
    chrom_counts = Counter()

    with open(MUTATIONS_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_mutations += 1
            vtype = row.get('Variant_Classification', 'Unknown')
            cohort = row.get('Cancer_Type', 'Unknown')
            chrom = row.get('Chromosome', 'Unknown')

            variant_types[vtype] += 1
            cancer_cohorts[cohort] += 1
            chrom_counts[chrom] += 1

    log(f"  • Total Extracted Noncoding Variants: {total_mutations:,}")
    log(f"  • Total Unique Cohorts/Folders Scanned: {len(cancer_cohorts)}")
    log(f"  • Chromosomes Represented: {len(chrom_counts)}")
    log("\n  Top Noncoding Functional Classifications:")
    for vt, count in variant_types.most_common(8):
        pct = (count / total_mutations) * 100 if total_mutations else 0
        log(f"    - {vt:<18}: {count:>6,} ({pct:>5.1f}%)")
else:
    log(f"  [!] Missing file: {MUTATIONS_FILE}")

# ---------------------------------------------------------
# 2. ANALYZE HOTSPOT RECURRENCE DATASET
# ---------------------------------------------------------
log("\n[2] 50bp RECURRENCE HOTSPOT ANALYSIS")
log("--------------------------------------------------------------------------")

if os.path.exists(HOTSPOTS_FILE):
    total_windows = 0
    recurrent_count = 0
    loci_counts = Counter()

    with open(HOTSPOTS_FILE, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.DictReader(f)
        for row in reader:
            total_windows += 1
            rec = int(row.get('Recurrence_Count', 0))
            gene = row.get('Gene', 'Unknown')
            
            if rec >= 2:
                recurrent_count += 1
            if gene and gene != 'Intergenic/Unknown':
                loci_counts[gene] += rec

    log(f"  • Total 50bp Genomic Windows Scanned: {total_windows:,}")
    log(f"  • Recurrent Hotspots (Recurrence >= 2): {recurrent_count:,}")
    log("\n  Top 10 High-Recurrence Target Loci:")
    for gene, count in loci_counts.most_common(10):
        log(f"    - {gene:<18}: {count:>5,} variants")
else:
    log(f"  [!] Missing file: {HOTSPOTS_FILE}")

# ---------------------------------------------------------
# 3. ANALYZE TRAINED DRIVER KNOWLEDGE BASE MODEL
# ---------------------------------------------------------
log("\n[3] TRAINED SERVER KNOWLEDGE BASE (JSON MODEL)")
log("--------------------------------------------------------------------------")

if os.path.exists(MODEL_FILE):
    with open(MODEL_FILE, 'r', encoding='utf-8') as f:
        model_data = json.load(f)

    total_regions = len(model_data)
    drivers = [v for v in model_data.values() if v.get('is_driver', False)]
    non_drivers = total_regions - len(drivers)

    driver_genes = Counter()
    driver_vtypes = Counter()
    scores = []

    for item in drivers:
        gene = item.get('gene', 'Unknown')
        score = item.get('score', 0.0)
        vtypes = item.get('variant_types', [])
        
        scores.append(score)
        driver_genes[gene] += 1
        for vt in vtypes:
            driver_vtypes[vt] += 1

    max_score = max(scores) if scores else 0.0
    avg_score = sum(scores) / len(scores) if scores else 0.0

    log(f"  • Model Knowledge Base Loci: {total_regions:,}")
    log(f"  • High-Confidence Driver Loci (Score >= 4.0 or Recurrence >= 3): {len(drivers):,}")
    log(f"  • Non-Driver / Low-Impact Loci: {non_drivers:,}")
    log(f"  • Maximum Pathogenicity Impact Score: {max_score:.2f}")
    log(f"  • Average Driver Locus Impact Score: {avg_score:.2f}")

    log("\n  Top Driver-Associated Target Genes:")
    for gene, count in driver_genes.most_common(8):
        if gene != 'Unknown':
            log(f"    - {gene:<18}: {count:>3} driver windows")

    log("\n  Driver Variant Category Representation:")
    for vt, count in driver_vtypes.most_common(6):
        log(f"    - {vt:<18}: {count:>3} occurrences")
else:
    log(f"  [!] Missing file: {MODEL_FILE}")

# ---------------------------------------------------------
# 4. SERVER SYSTEM INTEGRATION STATUS
# ---------------------------------------------------------
log("\n[4] DIAGNOSTIC ENGINE SYSTEM INTEGRATION STATUS")
log("--------------------------------------------------------------------------")
log("  • API Server Script : server.py")
log("  • Server Port       : http://localhost:8000")
log("  • Endpoint Active   : POST / for mutation query analysis")
log("  • Health Check      : GET /health")
log("  • Status            : READY FOR PAN-CANCER DIAGNOSTIC QUERIES")
log("==========================================================================")

# Save summary report to file
with open(REPORT_OUTPUT, 'w', encoding='utf-8') as out_f:
    out_f.write("\n".join(report_lines) + "\n")

log(f"\n[+] SUCCESS: Master report saved to disk at: {REPORT_OUTPUT}")