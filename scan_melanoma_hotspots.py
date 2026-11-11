"""
scan_melanoma_hotspots.py
Optimized Genome-Wide Noncoding Hotspot Discovery Engine for Melanoma.
Uses 100kb genomic binning to process 220,000+ hotspots in seconds.
"""

from pathlib import Path
import pandas as pd
import math
import time

RESULTS_DIR = Path(r"D:\noncoding mutation finding surver\results")
MELANOMA_FILE = RESULTS_DIR / "pcawg_melanoma_wgs_snvs.csv.gz"
OUTPUT_CSV = RESULTS_DIR / "tables" / "melanoma_noncoding_hotspots_ranked.csv"

def compute_poisson(k: int, lambd: float) -> float:
    if lambd <= 0:
        return 1.0 if k == 0 else 0.0
    prob_less = sum((math.exp(-lambd) * (lambd ** i)) / math.factorial(i) for i in range(k))
    return max(0.0, 1.0 - prob_less)

if not MELANOMA_FILE.exists():
    print(f"Error: {MELANOMA_FILE.name} not found.")
    exit(1)

print(f"Loading Melanoma WGS dataset ({MELANOMA_FILE.name})...")
start_time = time.time()
df = pd.read_csv(MELANOMA_FILE, low_memory=False)

# Clean up chromosomes and positions for fast grouping
df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False).str.strip()
df["Start_Position"] = pd.to_numeric(df["Start_Position"], errors="coerce").fillna(0).astype(int)

# 1. Identify Recurrent Hotspots
hotspots = df.groupby(["Chromosome", "Start_Position"])["patient"].nunique().reset_index(name="Mutated_Donors")
recurrent = hotspots[hotspots["Mutated_Donors"] >= 2].copy()
print(f"Identified {len(recurrent):,} recurrent somatic positions (N >= 2 donors).")

# Filter known coding exclusions
coding_exclusions = [("7", 140453136)] # BRAF V600E
for chrom, pos in coding_exclusions:
    recurrent = recurrent[~((recurrent["Chromosome"] == chrom) & (recurrent["Start_Position"] == pos))]

# 2. OPTIMIZATION: Pre-calculate background mutation rates in 100kb bins
print("Pre-calculating regional background rates...")
df["bin_100k"] = df["Start_Position"] // 100000
bg_counts = df.groupby(["Chromosome", "bin_100k"]).size().to_dict()

results = []
total_donors = 70

print("Evaluating Poisson significance for all hotspots (this will now take seconds)...")
for _, row in recurrent.iterrows():
    chrom = str(row["Chromosome"])
    pos = int(row["Start_Position"])
    n_mut = int(row["Mutated_Donors"])

    # Fast O(1) dictionary lookup for background mutations in that specific 100kb window
    bin_idx = pos // 100000
    regional_muts = bg_counts.get((chrom, bin_idx), 0)
    
    lambd = (regional_muts / 100000.0) * 51 if regional_muts > 0 else 0.05
    pval = compute_poisson(n_mut, lambd)
    fold_enrich = round(n_mut / max(lambd, 0.001), 2)

    results.append({
        "Chromosome": chrom,
        "Start_Position": pos,
        "Mutated_Donors": n_mut,
        "Cohort_Frequency_%": round((n_mut / total_donors) * 100, 2),
        "Expected_Baseline": round(lambd, 4),
        "Fold_Enrichment": fold_enrich,
        "Poisson_P_Value": pval
    })

res_df = pd.DataFrame(results).sort_values(by=["Mutated_Donors", "Fold_Enrichment"], ascending=[False, False])
res_df["Poisson_P_Value"] = res_df["Poisson_P_Value"].apply(lambda x: f"{x:.2e}")

OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
res_df.to_csv(OUTPUT_CSV, index=False)

print("\n" + "=" * 70)
print(" TOP NONCODING HOTSPOT CANDIDATES IN MELANOMA")
print("=" * 70)
print(res_df.head(10).to_string(index=False))
print("=" * 70)
print(f"\nCompleted in {round(time.time() - start_time, 1)} seconds.")
print(f"Saved full ranked candidate table to:\n  {OUTPUT_CSV}")