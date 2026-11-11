"""
03_explore_coding.py   (exploration only, not part of the noncoding project)

Which kinds of coding mutations occur in TCGA-SKCM melanoma exomes, and where?
Reads the MAF files directly, keeps one tumor sample per patient, then:
  1. counts coding mutations by class (missense, nonsense, silent, ...)
  2. shows the substitution spectrum (C>T, C>A, ...), which reveals UV damage
  3. lists the most recurrent protein changes (e.g. BRAF p.V600E)
  4. lists the most frequently mutated genes
  5. draws a lollipop plot (mutation position along the protein) for one gene

Outputs go to results/coding_explore/ inside your project folder.
"""
import glob
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # save figures to files; no pop-up windows
import matplotlib.pyplot as plt
import pandas as pd

# ========================= CHANGE THESE IF NEEDED ===========================
BASE       = Path(r"D:\noncoding mutation finding surver")  # <<< project folder
MAF_FOLDER = "maf_files"                                     # <<< folder created by 01_download.py
GENES      = ["BRAF", "NRAS", "TP53"]                        # <<< genes to draw lollipop plots for (any gene symbols)
# ============================================================================

out_dir = BASE / "results" / "coding_explore"
out_dir.mkdir(parents=True, exist_ok=True)
log_lines = []


def say(text=""):
    print(text)
    log_lines.append(str(text))


KEEP = {
    "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele",
    "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type",
    "Tumor_Sample_Barcode", "HGVSp_Short", "Protein_position",
}
CODING = {
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent",
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Translation_Start_Site",
}
SAMPLE_PRIORITY = {"01": 0, "06": 1, "03": 2, "02": 3, "05": 4}  # lower = preferred


def n_header_lines(path):
    n = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                n += 1
            else:
                break
    return n


# ------------------------------ read and filter -----------------------------
files = sorted(glob.glob(str(BASE / MAF_FOLDER / "**" / "*.maf.gz"), recursive=True))
say(f"MAF files found: {len(files)}")
if not files:
    raise SystemExit("No .maf.gz files found. Check BASE and MAF_FOLDER at the top of the script.")

frames = []
for i, f in enumerate(files, 1):
    try:
        df = pd.read_csv(f, sep="\t", skiprows=n_header_lines(f), compression="gzip",
                         usecols=lambda c: c in KEEP, dtype=str, low_memory=False)
        frames.append(df[df["Variant_Classification"].isin(CODING)])
    except Exception as e:
        say(f"  could not read {Path(f).name}: {str(e)[:80]}")
    if i % 100 == 0:
        print(f"  read {i}/{len(files)} files")

maf = pd.concat(frames, ignore_index=True)
for col in ("HGVSp_Short", "Protein_position"):
    if col not in maf.columns:
        say(f"Note: column {col} is not in these files, so sections that use it are skipped.")
        maf[col] = pd.NA

# one tumor sample per patient (same rule as script 02)
maf["patient"] = maf["Tumor_Sample_Barcode"].str[:12]
maf["sample_code"] = maf["Tumor_Sample_Barcode"].str[13:15]
samples = maf[["patient", "Tumor_Sample_Barcode", "sample_code"]].drop_duplicates()
samples["rank"] = samples["sample_code"].map(SAMPLE_PRIORITY).fillna(9)
chosen = (samples.sort_values(["patient", "rank", "Tumor_Sample_Barcode"])
                 .drop_duplicates("patient")[["Tumor_Sample_Barcode"]])
maf = maf.merge(chosen, on="Tumor_Sample_Barcode")
maf = maf.drop_duplicates(["Tumor_Sample_Barcode", "Chromosome", "Start_Position",
                           "Reference_Allele", "Tumor_Seq_Allele2"]).reset_index(drop=True)
n_patients = maf["patient"].nunique()
say(f"Patients: {n_patients}   Coding mutations: {len(maf):,}")

nonsilent = maf[maf["Variant_Classification"] != "Silent"]

# ------------------------- 1. counts by mutation class ----------------------
cls = maf["Variant_Classification"].value_counts()
say("\n1. Coding mutations by class:")
say(pd.DataFrame({"mutations": cls, "percent": (100 * cls / cls.sum()).round(2)}).to_string())

# ----------------------- 2. substitution spectrum ---------------------------
snv = maf[(maf["Variant_Type"] == "SNP")
          & maf["Reference_Allele"].isin(list("ACGT"))
          & maf["Tumor_Seq_Allele2"].isin(list("ACGT"))]
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
flip = snv["Reference_Allele"].isin(["A", "G"])   # report every change from the C or T strand
ref = snv["Reference_Allele"].where(~flip, snv["Reference_Allele"].map(comp))
alt = snv["Tumor_Seq_Allele2"].where(~flip, snv["Tumor_Seq_Allele2"].map(comp))
spectrum = (ref + ">" + alt).value_counts().reindex(
    ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]).fillna(0).astype(int)
say("\n2. Substitution spectrum of coding point mutations (C>T = UV signature in melanoma):")
say(pd.DataFrame({"count": spectrum, "percent": (100 * spectrum / spectrum.sum()).round(1)}).to_string())

# ------------------- 3. most recurrent protein changes ----------------------
prot = nonsilent.dropna(subset=["HGVSp_Short"])
prot = prot[prot["HGVSp_Short"].str.len() > 0]
recur = (prot.groupby(["Hugo_Symbol", "HGVSp_Short"])["patient"].nunique()
             .sort_values(ascending=False).head(25))
say("\n3. Most recurrent protein changes (patients carrying it, of "
    f"{n_patients}):")
say(pd.DataFrame({"patients": recur, "percent": (100 * recur / n_patients).round(1)}).to_string())

# ----------------------- 4. most frequently mutated genes -------------------
genes = nonsilent.groupby("Hugo_Symbol")["patient"].nunique().sort_values(ascending=False).head(25)
say("\n4. Genes with non-silent coding mutations in the most patients:")
say(pd.DataFrame({"patients": genes, "percent": (100 * genes / n_patients).round(1)}).to_string())
say("CAUTION: very large genes (TTN, MUC16, ...) top such lists because they offer more "
    "positions to mutate, and melanoma has a very high mutation rate. A frequently mutated "
    "gene is not automatically a cancer driver.")

# ---------------------------- 5. lollipop plots -----------------------------
COLORS = {
    "Missense_Mutation": "tab:green", "Nonsense_Mutation": "tab:red",
    "Frame_Shift_Del": "tab:purple", "Frame_Shift_Ins": "tab:purple",
    "In_Frame_Del": "tab:orange", "In_Frame_Ins": "tab:orange",
    "Silent": "lightgray", "Nonstop_Mutation": "tab:brown",
    "Translation_Start_Site": "tab:pink",
}
for gene in GENES:
    sub = maf[maf["Hugo_Symbol"] == gene].copy()
    sub["pos"] = pd.to_numeric(sub["Protein_position"].str.extract(r"(\d+)")[0], errors="coerce")
    sub = sub.dropna(subset=["pos"])
    if sub.empty:
        say(f"\n5. {gene}: no coding mutations with a protein position found, skipped.")
        continue
    length = pd.to_numeric(sub["Protein_position"].str.extract(r"/(\d+)")[0], errors="coerce").max()
    g = sub.groupby(["pos", "Variant_Classification"])["patient"].nunique().reset_index(name="patients")
    n_gene = sub[sub["Variant_Classification"] != "Silent"]["patient"].nunique()
    say(f"\n5. {gene}: {n_gene} patients ({100 * n_gene / n_patients:.1f}%) have a non-silent coding mutation. "
        "Top positions:")
    say(g.sort_values("patients", ascending=False).head(8).to_string(index=False))

    fig, ax = plt.subplots(figsize=(11, 3.8))
    for c, part in g.groupby("Variant_Classification"):
        ax.vlines(part["pos"], 0, part["patients"], color=COLORS.get(c, "black"), lw=1)
        ax.scatter(part["pos"], part["patients"], s=28, color=COLORS.get(c, "black"), label=c, zorder=3)
    if pd.notna(length):
        ax.set_xlim(0, length)
    ax.set_xlabel("Position along the protein (amino acid number)")
    ax.set_ylabel("Patients")
    ax.set_title(f"{gene}: coding mutations in TCGA-SKCM")
    ax.legend(fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(out_dir / f"lollipop_{gene}.png", dpi=200)
    plt.close(fig)

# ------------------------------ plots and files -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
cls.sort_values().plot(kind="barh", ax=axes[0])
axes[0].set_title("Coding mutations by class")
spectrum.plot(kind="bar", ax=axes[1], rot=0)
axes[1].set_title("Substitution spectrum (coding point mutations)")
plt.tight_layout()
plt.savefig(out_dir / "class_and_spectrum.png", dpi=200)
plt.close(fig)

maf.to_csv(out_dir / "coding_mutations.csv.gz", index=False)
(out_dir / "summary_coding.txt").write_text("\n".join(log_lines), encoding="utf-8")
say(f"\nDone. Outputs are in: {out_dir}")