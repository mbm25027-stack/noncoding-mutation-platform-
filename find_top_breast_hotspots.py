"""
find_top_breast_hotspots.py
Scans pcawg_breast_wgs_snvs.csv.gz to identify top recurrent noncoding coordinates.
"""

from pathlib import Path
import pandas as pd

results_dir = Path(r"D:\noncoding mutation finding surver\results")
breast_file = results_dir / "pcawg_breast_wgs_snvs.csv.gz"

if not breast_file.exists():
    print(f"File not found: {breast_file}")
else:
    print("Indexing Breast Cancer WGS SNVs to find top recurrent loci...")
    df = pd.read_csv(breast_file, low_memory=False)

    chrom_col = "Chromosome" if "Chromosome" in df.columns else df.columns[0]
    pos_col = "Start_Position" if "Start_Position" in df.columns else df.columns[1]
    patient_col = "patient" if "patient" in df.columns else ("Donor_ID" if "Donor_ID" in df.columns else df.columns[2])

    df["Chromosome"] = df[chrom_col].astype(str).str.replace("chr", "", case=False).str.strip()
    df["Start_Position"] = pd.to_numeric(df[pos_col], errors="coerce").fillna(0).astype(int)

    # Group by Chromosome and Position to find most mutated sites
    hotspots = (
        df.groupby(["Chromosome", "Start_Position"])[patient_col]
        .nunique()
        .reset_index(name="Mutated_Donors")
        .sort_values(by="Mutated_Donors", ascending=False)
    )

    print("\n" + "=" * 60)
    print(" TOP 10 RECURRENT SOMATIC LOCI IN BREAST CANCER COHORT")
    print("=" * 60)
    print(hotspots.head(10).to_string(index=False))
    print("=" * 60)