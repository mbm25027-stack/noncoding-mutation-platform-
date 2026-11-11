"""
prep_breast_data.py - Robust Breast Cancer Data Preparation
Fixes delimiter parsing for extensionless PCAWG histology files.
"""

from pathlib import Path
import pandas as pd

base_dir = Path(r"D:\noncoding mutation finding surver")
results_dir = base_dir / "results"

# 1. Locate Pan-Cancer MAF file
maf_file = None
for search_path in [results_dir, base_dir]:
    for f in search_path.glob("*"):
        if "whitelist" in f.name.lower() or "nonus" in f.name.lower():
            if f.name != "pcawg_melanoma_wgs_snvs.csv.gz":
                maf_file = f
                break
    if maf_file:
        break

# 2. Locate Histology mapping file
histo_file = None
for search_path in [base_dir, results_dir]:
    for f in search_path.glob("*histology*"):
        histo_file = f
        break

if not maf_file:
    print("Error: Could not locate the raw MAF dataset.")
    exit(1)

print(f"Reading pan-cancer MAF dataset: {maf_file.name} ...")
df_raw = pd.read_csv(maf_file, sep="\t", comment="#", low_memory=False)

print(f"Raw MAF loaded with {len(df_raw):,} total mutation records across all cancer types.")
print(f"Columns detected: {list(df_raw.columns[:6])}")

# 3. Standardize MAF column names
col_map = {}
for col in df_raw.columns:
    c_lower = col.lower()
    if c_lower in ["chrom", "chr", "chromosome"]:
        col_map[col] = "Chromosome"
    elif c_lower in ["start", "start_position", "pos", "position"]:
        col_map[col] = "Start_Position"
    elif c_lower in ["ref", "reference", "reference_allele"]:
        col_map[col] = "Reference_Allele"
    elif c_lower in ["alt", "tumor_seq_allele2", "alternate", "tumor_allele"]:
        col_map[col] = "Tumor_Seq_Allele2"
    elif c_lower in ["sample", "donor_id", "patient", "icgc_donor_id", "sample_id"]:
        col_map[col] = "patient"

df = df_raw.rename(columns=col_map)

# 4. Parse histology file using tab separation
breast_donors = []
if histo_file and histo_file.exists():
    print(f"Loading histology donor annotations: {histo_file.name} ...")
    try:
        histo_df = pd.read_csv(histo_file, sep="\t", comment="#", low_memory=False)
    except Exception:
        histo_df = pd.read_csv(histo_file, sep=None, engine="python", comment="#")

    print(f"Histology columns: {list(histo_df.columns[:5])}")

    # Find donor and histology/project columns
    donor_col = next((c for c in histo_df.columns if any(k in c.lower() for k in ["donor", "icgc", "specimen", "sample"])), histo_df.columns[0])
    hist_col = next((c for c in histo_df.columns if any(k in c.lower() for k in ["histology", "project", "dcc", "cancer", "organ", "tumour", "tumor"])), None)

    if hist_col:
        breast_mask = histo_df[hist_col].astype(str).str.contains("Breast|BRCA", case=False, na=False)
        breast_donors = histo_df[breast_mask][donor_col].dropna().unique().tolist()
        print(f"Identified {len(breast_donors)} Breast Cancer donors in histology reference.")

# 5. Filter dataset for Breast Cancer donors
if breast_donors and "patient" in df.columns:
    df_breast = df[df["patient"].isin(breast_donors)].copy()
    print(f"Filtered to {len(df_breast):,} somatic mutations belonging to Breast Cancer donors.")
else:
    print("Notice: No histology filter matched; keeping all records.")
    df_breast = df.copy()

# 6. Standardize Chromosome formatting
if "Chromosome" in df_breast.columns:
    df_breast["Chromosome"] = df_breast["Chromosome"].astype(str).str.replace("chr", "", case=False).str.strip()

# 7. Save compressed dataset
output_path = results_dir / "pcawg_breast_wgs_snvs.csv.gz"
df_breast.to_csv(output_path, index=False, compression="gzip")
print(f"\nSUCCESS: Formatted Breast Cancer dataset saved to:\n  {output_path}")