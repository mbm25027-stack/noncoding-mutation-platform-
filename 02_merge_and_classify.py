"""
02_merge_and_classify.py

Merge the TCGA-SKCM exome MAF files downloaded by 01_download.py, keep one
tumor sample per patient, label every mutation by region class, and report how
much noncoding signal the exome data holds.

Outputs (written to the results folder):
  all_mutations.csv.gz         every mutation, one sample per patient
  noncoding_mutations.csv.gz   only the noncoding classes
  region_counts.csv            counts per Variant_Classification
  region_counts.png            bar charts of the same
  summary.txt                  everything printed to the screen (for your README)
"""
import glob
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # save figures to a file instead of opening a window
import matplotlib.pyplot as plt
import pandas as pd

# ========================= CHANGE THESE IF NEEDED ===========================
BASE       = Path(r"D:\noncoding mutation finding surver")  # <<< project folder
MAF_FOLDER = "maf_files"   # <<< folder created by 01_download.py
OUT_FOLDER = "results"     # <<< where the outputs are written
# ============================================================================

out_dir = BASE / OUT_FOLDER
out_dir.mkdir(exist_ok=True)
log_lines = []


def say(text=""):
    print(text)
    log_lines.append(str(text))


# Columns we need from each MAF (the files have ~120 columns, we load only these)
KEEP = {
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
    "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification",
    "Variant_Type", "Tumor_Sample_Barcode", "NCBI_Build", "HGVSc",
    "t_depth", "t_alt_count",
}

# Region class of every Variant_Classification value found in GDC MAF files
REGION_GROUP = {
    "Missense_Mutation": "coding", "Nonsense_Mutation": "coding",
    "Nonstop_Mutation": "coding", "Silent": "coding",
    "Frame_Shift_Del": "coding", "Frame_Shift_Ins": "coding",
    "In_Frame_Del": "coding", "In_Frame_Ins": "coding",
    "Translation_Start_Site": "coding",
    "3'UTR": "UTR", "5'UTR": "UTR",
    # Splice sites/regions lie next to exons: counted as noncoding here, although
    # they can still change the protein. Keep this in mind when you interpret them.
    "Splice_Site": "splice", "Splice_Region": "splice",
    "Intron": "intron",
    "3'Flank": "flank", "5'Flank": "flank",
    "RNA": "ncRNA",
    "IGR": "intergenic",
}
NONCODING_GROUPS = {"UTR", "splice", "intron", "flank", "ncRNA", "intergenic"}

# TCGA sample type codes (characters 14-15 of the barcode): 01 primary tumor,
# 06 metastatic, 03 blood-derived (leukemia), 02 recurrent, 05 additional new primary
SAMPLE_PRIORITY = {"01": 0, "06": 1, "03": 2, "02": 3, "05": 4}  # lower = preferred


def n_header_lines(path):
    """Count the leading '#' comment lines (GDC MAFs start with '#version ...')."""
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                n += 1
            else:
                break
    return n


# ------------------------------ 1. read files ------------------------------
files = sorted(glob.glob(str(BASE / MAF_FOLDER / "**" / "*.maf.gz"), recursive=True))
say(f"MAF files found: {len(files)}")
if not files:
    raise SystemExit("No .maf.gz files found. Check BASE and MAF_FOLDER at the top of the script.")

frames, failed = [], []
for i, f in enumerate(files, 1):
    try:
        df = pd.read_csv(
            f, sep="\t", skiprows=n_header_lines(f), compression="gzip",
            usecols=lambda c: c in KEEP, dtype=str, low_memory=False,
        )
        frames.append(df)
    except Exception as e:  # keep going if one file is damaged
        failed.append((Path(f).name, str(e)[:80]))
    if i % 100 == 0:
        print(f"  read {i}/{len(files)} files")

if failed:
    say(f"WARNING: {len(failed)} file(s) could not be read:")
    for name, err in failed[:10]:
        say(f"  {name}: {err}")

maf = pd.concat(frames, ignore_index=True)
say(f"Rows read: {len(maf):,}")

missing_cols = KEEP - set(maf.columns)
if missing_cols:
    say(f"Note: these optional columns are not in the files: {sorted(missing_cols)}")

maf["Start_Position"] = pd.to_numeric(maf["Start_Position"], errors="coerce")

# ------------------ 2. one tumor sample per patient ------------------------
maf["patient"] = maf["Tumor_Sample_Barcode"].str[:12]
maf["sample_code"] = maf["Tumor_Sample_Barcode"].str[13:15]

samples = maf[["patient", "Tumor_Sample_Barcode", "sample_code"]].drop_duplicates()
samples["rank"] = samples["sample_code"].map(SAMPLE_PRIORITY).fillna(9)
n_multi = int((samples.groupby("patient").size() > 1).sum())

chosen = (samples.sort_values(["patient", "rank", "Tumor_Sample_Barcode"])
                 .drop_duplicates("patient")[["Tumor_Sample_Barcode"]])
maf = maf.merge(chosen, on="Tumor_Sample_Barcode")

key = ["Tumor_Sample_Barcode", "Chromosome", "Start_Position",
       "Reference_Allele", "Tumor_Seq_Allele2"]
maf = maf.drop_duplicates(key).reset_index(drop=True)

say(f"Patients: {maf['patient'].nunique()}  "
    f"(patients with more than one sample in the files: {n_multi}; "
    f"kept one each, preferring primary '01', then metastatic '06')")
say("Sample types kept:")
say(maf.drop_duplicates("patient")["sample_code"].value_counts().to_string())
say(f"Mutations after de-duplication: {len(maf):,}")

# ------------------------------ 3. sanity checks ---------------------------
say("\nGenome build(s) in the files:")
say(maf["NCBI_Build"].value_counts(dropna=False).to_string())
say("\nVariant types:")
say(maf["Variant_Type"].value_counts(dropna=False).to_string())

# ----------------------------- 4. region classes ---------------------------
maf["region_group"] = maf["Variant_Classification"].map(REGION_GROUP).fillna("other")
maf["is_noncoding"] = maf["region_group"].isin(NONCODING_GROUPS)

counts = maf["Variant_Classification"].value_counts()
table = counts.rename("mutations").to_frame()
table["percent"] = (100 * table["mutations"] / table["mutations"].sum()).round(2)
table["region_group"] = table.index.map(lambda c: REGION_GROUP.get(c, "other"))
table.to_csv(out_dir / "region_counts.csv")

say("\nMutations by Variant_Classification:")
say(table.to_string())

group_counts = maf["region_group"].value_counts()
say("\nMutations by region group:")
say(group_counts.to_string())
say(f"\nNoncoding fraction of all mutations: {maf['is_noncoding'].mean():.3f}")

unknown = sorted(set(maf["Variant_Classification"].dropna()) - set(REGION_GROUP))
if unknown:
    say(f"Classes not in the mapping (labelled 'other'): {unknown}")

# --------------------- 5. noncoding RNAs (your focus) ----------------------
rna = maf[maf["Variant_Classification"] == "RNA"]
say(f"\nMutations in noncoding RNA genes ('RNA' class): {len(rna):,}")
if len(rna):
    say("Most frequently mutated ncRNA genes (mutations):")
    say(rna["Hugo_Symbol"].value_counts().head(15).to_string())

# ---------- 6. preview of recurrent noncoding positions (NOT hotspots) -----
nc = maf[maf["is_noncoding"]]
recur = (nc.groupby(["Chromosome", "Start_Position", "Reference_Allele",
                     "Tumor_Seq_Allele2", "Variant_Classification", "Hugo_Symbol"])["patient"]
           .nunique().sort_values(ascending=False))
recur = recur[recur >= 3]
say("\nNoncoding positions mutated in >= 3 patients (top 20):")
say(recur.head(20).to_string() if len(recur) else "  none")
say("CAUTION: this is only a preview. Melanoma has a very high UV mutation rate, so many "
    "positions recur by chance. Real hotspot calling needs a background mutation model.")

# ------------------------------- 7. figure ---------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
counts.sort_values().plot(kind="barh", ax=axes[0])
axes[0].set_title("TCGA-SKCM exome: mutations by class")
axes[0].set_xlabel("Number of mutations")
group_counts.plot(kind="bar", ax=axes[1], rot=30)
axes[1].set_title("Mutations by region group")
axes[1].set_ylabel("Number of mutations")
plt.tight_layout()
plt.savefig(out_dir / "region_counts.png", dpi=200)
plt.close(fig)

# ------------------------------ 8. save tables -----------------------------
maf.to_csv(out_dir / "all_mutations.csv.gz", index=False)
nc.to_csv(out_dir / "noncoding_mutations.csv.gz", index=False)
(out_dir / "summary.txt").write_text("\n".join(log_lines), encoding="utf-8")

say(f"\nDone. Outputs are in: {out_dir}")