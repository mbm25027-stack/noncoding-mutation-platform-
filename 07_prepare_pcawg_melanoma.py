"""
07_prepare_pcawg_melanoma.py

Parse open-access PCAWG whole-genome somatic mutations downloaded from UCSC Xena:
  1. Identifies melanoma donors from tumour histology metadata.
  2. Extracts all whole-genome somatic SNVs for melanoma patients.
  3. Verifies TERT promoter positive control hotspots (C228T and C250T).

Outputs:
  - results/pcawg_melanoma_wgs_snvs.csv.gz
"""

from pathlib import Path
import pandas as pd

# ========================= PATH CONFIGURATION ===========================
BASE = Path(r"D:\noncoding mutation finding surver")
OUT_DIR = BASE / "results"
OUT_DIR.mkdir(exist_ok=True)
# ========================================================================


def find_pcawg_file(patterns, label):
    for pat in patterns:
        hits = list(BASE.rglob(pat))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"Could not locate {label} using patterns: {patterns}")


# 1. Locate downloaded files
mut_path = find_pcawg_file(
    ["*whitelist*.nonUS*", "*whitelist*.maf*", "*.snv_mnv_indel*"],
    "PCAWG WGS Mutation File",
)
hist_path = find_pcawg_file(
    ["*histology*.tsv*", "*specimen_histology*"],
    "PCAWG Histology Metadata File",
)

print(f"Mutation File:  {mut_path.name}")
print(f"Histology File: {hist_path.name}\n")

# 2. Identify Melanoma Donors
hist = pd.read_csv(hist_path, sep="\t", low_memory=False)

donor_col = next(c for c in hist.columns if any(k in c.lower() for k in ["donor", "sample", "icgc"]))
project_col = next((c for c in hist.columns if any(k in c.lower() for k in ["project", "histology", "tumour", "abbreviation"])), None)

print(f"Donor Column: {donor_col} | Project/Histology Column: {project_col}")

melanoma_mask = hist[project_col].astype(str).str.contains("MELA|Melanoma|Skin", case=False, na=False)
melanoma_donors = set(hist.loc[melanoma_mask, donor_col].dropna().unique())

print(f"Total melanoma donors found in metadata: {len(melanoma_donors)}")

# 3. Stream & Filter Whole-Genome Mutation File
print("\nReading PCAWG whole-genome mutations (filtering for melanoma donors)...")
chunk_size = 500000
melanoma_chunks = []

for chunk in pd.read_csv(mut_path, sep="\t", chunksize=chunk_size, low_memory=False):
    sample_col = next(c for c in chunk.columns if any(k in c.lower() for k in ["sample", "donor", "patient"]))
    sub = chunk[chunk[sample_col].isin(melanoma_donors)].copy()
    if not sub.empty:
        melanoma_chunks.append(sub)

if not melanoma_chunks:
    raise SystemExit("No exact sample ID matches found between mutation file and histology donor list.")

melanoma_muts = pd.concat(melanoma_chunks, ignore_index=True)
print(f"Melanoma WGS mutations extracted: {len(melanoma_muts):,}")

# Standardize column names dynamically
col_map = {}
for c in melanoma_muts.columns:
    cl = c.lower().strip()
    if cl in ["chr", "chrom", "chromosome", "#chrom"] or ("chr" in cl and "start" not in cl and "end" not in cl):
        col_map[c] = "Chromosome"
    elif "start" in cl:
        col_map[c] = "Start_Position"
    elif "ref" in cl:
        col_map[c] = "Reference_Allele"
    elif "alt" in cl:
        col_map[c] = "Tumor_Seq_Allele2"
    elif any(k in cl for k in ["sample", "donor", "patient"]):
        col_map[c] = "patient"

melanoma_muts = melanoma_muts.rename(columns=col_map)
print(f"Mapped columns: {list(melanoma_muts.columns[:8])}")

# Clean Chromosome & Position formats
melanoma_muts["Chromosome"] = melanoma_muts["Chromosome"].astype(str).str.replace("chr", "")
melanoma_muts["Start_Position"] = pd.to_numeric(melanoma_muts["Start_Position"], errors="coerce")

# 4. TERT Promoter Positive Control Sanity Check (GRCh37 / hg19 coordinates)
print("\n--- SANITY CHECK: TERT Promoter Hotspots (chr5) ---")
pos = melanoma_muts["Start_Position"]
chrom_5 = melanoma_muts["Chromosome"] == "5"

tert_region = melanoma_muts[chrom_5 & pos.between(1295200, 1295260)]
c228t = melanoma_muts[chrom_5 & pos.between(1295227, 1295229)]
c250t = melanoma_muts[chrom_5 & pos.between(1295249, 1295251)]

print(f"Mutations near C228T (chr5:1295228): {len(c228t)} variants")
print(f"Mutations near C250T (chr5:1295250): {len(c250t)} variants")
print(f"Total TERT promoter region mutations found (1295200-1295260): {len(tert_region)}")

# 5. Save output
out_file = OUT_DIR / "pcawg_melanoma_wgs_snvs.csv.gz"
melanoma_muts.to_csv(out_file, index=False, compression="gzip")
print(f"\nSaved extracted melanoma WGS mutations to: {out_file}")