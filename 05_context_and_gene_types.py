"""
05_context_and_gene_types.py

Adds two things to every mutation and summarises them:

  1. Trinucleotide context. The reference genome supplies the base before and
     after each point mutation, which gives the 96-class mutation spectrum
     (for example A[C>T]G = a C>T change with A on the 5' side and G on the 3' side).
  2. GENCODE gene type. The mutated gene is classed as protein_coding, lncRNA,
     miRNA, other small ncRNA, pseudogene, and so on.

Why this step matters for the noncoding project
  * The 96-class spectrum is the mutational fingerprint. It is the core input
    of the background mutation model that hotspot calling needs (script 06).
  * Gene types replace the unreliable 'RNA' label with real lncRNA / miRNA
    annotation, so you can study true noncoding RNAs.
  * Matching the reference base to the MAF reference allele proves that the
    genome file and the coordinates agree, which protects every later step.

Inputs
    results/all_mutations.csv.gz            from 02_merge_and_classify.py
    GRCh38.d1.vd1.fa                        GDC reference genome (found automatically)
    gencode.v36.annotation.gtf              GENCODE v36 annotation (found automatically)
Needs:  pip install pyfaidx   (use the same Python that runs your scripts)

Outputs
    results/context/mutations_with_context.csv.gz
    results/tables/13_* ... 18_*  (CSV)      results/figures/fig10_*, fig11_*  (PNG)
"""
import gzip
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # write files, never open windows
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from pyfaidx import Fasta
except ImportError:
    raise SystemExit("pyfaidx is not installed. Run:  python -m pip install pyfaidx  "
                     "(with the same Python you use to run this script).")

# ========================= CHANGE THESE IF NEEDED ===========================
BASE = Path(r"D:\noncoding mutation finding surver")   # <<< project folder
FASTA = None   # <<< None = find GRCh38.d1.vd1.fa automatically, or give a full path in r"..."
GTF = None     # <<< None = find gencode.v36.annotation.gtf automatically, or give a full path
# ============================================================================

RES = BASE / "results"
TAB = RES / "tables"
FIG = RES / "figures"
CTX = RES / "context"

SPECTRUM = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BASES = "ACGT"
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
CHANNELS_96 = [f"{five}[{sub}]{three}" for sub in SPECTRUM for five in BASES for three in BASES]
SUB_COLORS = {"C>A": "#03BCEE", "C>G": "#010101", "C>T": "#E32926",
              "T>A": "#999999", "T>C": "#A1CE63", "T>G": "#EBC6C4"}
GROUP_ORDER = ["protein_coding", "lncRNA", "miRNA", "other_small_ncRNA", "pseudogene", "other", "unmatched"]
GROUP_COLORS = {"protein_coding": "#1F3864", "lncRNA": "#6A4C93", "miRNA": "#E76F51",
                "other_small_ncRNA": "#F4A261", "pseudogene": "#8D99AE", "other": "#CCCCCC",
                "unmatched": "#999999"}


# ------------------------------ small helpers -------------------------------
def save_table(df, name, index=True):
    path = TAB / f"{name}.csv"
    df.to_csv(path, index=index)
    print(f"  table   {path.name}")


def save_fig(fig, name):
    path = FIG / f"{name}.png"
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  figure  {path.name}")


def find_file(configured, names, what):
    """Return a real file. Windows extraction often creates a FOLDER with the
    file's name and the file inside it, so search recursively for files only."""
    if configured:
        p = Path(configured)
        if p.is_file():
            return p
        raise SystemExit(f"{what}: {p} is not a file.")
    for name in names:
        hits = [p for p in BASE.rglob(name) if p.is_file()]
        if hits:
            return sorted(hits, key=lambda p: -p.stat().st_size)[0]
    raise SystemExit(f"Could not find {what} under {BASE}. Looked for: {names}. "
                     "Set FASTA / GTF at the top of the script.")


# -------------------------- GENCODE gene types ------------------------------
def biotype_group(t):
    """Collapse the ~40 GENCODE gene types into a few groups."""
    if pd.isna(t):
        return "unmatched"
    if "pseudogene" in t:
        return "pseudogene"
    if t == "protein_coding" or (t.startswith(("IG_", "TR_")) and t.endswith("_gene")):
        return "protein_coding"
    if t == "lncRNA":
        return "lncRNA"
    if t == "miRNA":
        return "miRNA"
    if t in {"snRNA", "snoRNA", "scaRNA", "misc_RNA", "rRNA", "Mt_rRNA", "Mt_tRNA",
             "sRNA", "scRNA", "vault_RNA", "ribozyme"}:
        return "other_small_ncRNA"
    return "other"


def load_gene_table(gtf_path):
    """Read only the 'gene' lines of the GENCODE GTF (about 60,000 lines).

    WHY: gives each gene symbol its official gene type, so noncoding RNAs can be
    defined from GENCODE rather than from the MAF's 'RNA' label."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    pat = {k: re.compile(rf'{k} "([^"]+)"') for k in ("gene_id", "gene_name", "gene_type", "gene_biotype")}
    rows = []
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "\tgene\t" not in line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            attr = f[8]
            gid = pat["gene_id"].search(attr)
            name = pat["gene_name"].search(attr)
            typ = pat["gene_type"].search(attr) or pat["gene_biotype"].search(attr)
            rows.append((gid.group(1) if gid else None, name.group(1) if name else None,
                         typ.group(1) if typ else None))
    genes = pd.DataFrame(rows, columns=["gene_id", "gene_name", "gene_type"]).dropna(subset=["gene_name"])
    genes["gene_group"] = genes["gene_type"].map(biotype_group)
    # keep one row per symbol, preferring ordinary (non-PAR) gene ids
    genes["is_par"] = genes["gene_id"].fillna("").str.contains("_PAR_")
    genes = genes.sort_values("is_par").drop_duplicates("gene_name")
    print(f"  GENCODE genes read: {len(genes):,}")
    return genes.set_index("gene_name")[["gene_id", "gene_type", "gene_group"]]


# ---------------------------- reference context -----------------------------
def add_context(df, genome):
    """Look up the reference base before, at and after each point mutation.

    WHY: the surrounding bases decide how likely a change is (UV damage hits
    dipyrimidines, ageing hits CpG). Without context a background model would
    mistake predictable mutation-prone sites for hotspots."""
    df = df.copy()
    df["trinuc"] = ""
    chroms_in_genome = set(genome.keys())
    missing = sorted(set(df["Chromosome"].astype(str)) - chroms_in_genome)
    if missing:
        print(f"  note: {len(missing)} chromosome name(s) in the mutations are not in the FASTA "
              f"(first few: {missing[:5]}); those rows get no context")

    snv_mask = (df["Variant_Type"] == "SNP") & df["Reference_Allele"].isin(list(BASES)) \
        & df["Tumor_Seq_Allele2"].isin(list(BASES))
    work = df[snv_mask]
    for chrom, idx in work.groupby("Chromosome").groups.items():
        if chrom not in chroms_in_genome:
            continue
        raw = genome[chrom][:].seq.encode("ascii")           # whole chromosome, loaded once
        seq = np.frombuffer(raw, dtype=np.uint8) & 0xDF      # force upper case
        n = len(seq)
        pos = df.loc[idx, "Start_Position"].to_numpy(dtype=np.int64)  # 1-based
        ok = (pos >= 2) & (pos <= n - 1)
        pos = np.where(ok, pos, 2)
        tri = np.ascontiguousarray(np.stack([seq[pos - 2], seq[pos - 1], seq[pos]], axis=1))
        tri = tri.view("S3")[:, 0].astype("U3")
        tri = np.where(ok, tri, "")
        df.loc[idx, "trinuc"] = tri
        del raw, seq
        print(f"    context added: {chrom}", end="\r")
    print(" " * 40, end="\r")

    valid = df["trinuc"].str.fullmatch(r"[ACGT]{3}", na=False)
    df["ref_match"] = valid & (df["trinuc"].str[1] == df["Reference_Allele"])
    df["has_context"] = valid & snv_mask
    return df


def fold_label(ref, alt, tri):
    """Express a change on the pyrimidine (C or T) strand: 'A[C>T]G'."""
    if ref in "CT":
        return f"{tri[0]}[{ref}>{alt}]{tri[2]}"
    return f"{COMP[tri[2]]}[{COMP[ref]}>{COMP[alt]}]{COMP[tri[0]]}"


def add_spectrum_labels(df):
    use = df["has_context"] & df["ref_match"]
    labels = pd.Series("", index=df.index)
    labels[use] = [fold_label(r, a, t) for r, a, t in
                   zip(df.loc[use, "Reference_Allele"], df.loc[use, "Tumor_Seq_Allele2"], df.loc[use, "trinuc"])]
    df["sub96"] = labels
    df["sub6"] = df["sub96"].str[2:5]
    return df


# ---------------------------------- tables ----------------------------------
def table_reference_check(df):
    """Table: does the reference genome agree with the MAF reference alleles?
    WHY: if the match rate is not ~100%, the genome file or its coordinates are
    wrong, and every context-based result would be meaningless."""
    n_snv = int((df["Variant_Type"] == "SNP").sum())
    n_ctx = int(df["has_context"].sum())
    n_match = int(df["ref_match"].sum())
    t = pd.DataFrame({"value": [n_snv, n_ctx, n_match, round(100 * n_match / max(n_ctx, 1), 3)]},
                     index=["point mutations", "with reference context", "reference base matches MAF",
                            "match rate (%)"])
    save_table(t, "13_reference_check")
    rate = 100 * n_match / max(n_ctx, 1)
    print(f"  reference match rate: {rate:.3f}%")
    if rate < 99:
        print("  WARNING: match rate below 99%. Check that the FASTA is GRCh38.d1.vd1 and that "
              "chromosome names agree. Do not trust the context results yet.")


def spectrum96(df):
    use = df[df["sub96"] != ""]
    cod = use[use["region_group"] == "coding"]["sub96"].value_counts().reindex(CHANNELS_96).fillna(0).astype(int)
    nc = use[use["is_noncoding"]]["sub96"].value_counts().reindex(CHANNELS_96).fillna(0).astype(int)
    return cod, nc


def table_spectrum96(cod, nc):
    """Table: the 96-class mutation spectrum, coding and noncoding.
    WHY: this is the mutational fingerprint your background model will use."""
    t = pd.DataFrame({"coding_count": cod, "coding_percent": (100 * cod / max(cod.sum(), 1)).round(3),
                      "noncoding_count": nc, "noncoding_percent": (100 * nc / max(nc.sum(), 1)).round(3)})
    save_table(t, "14_spectrum96")


def table_uv_check(cod, nc):
    """Table: rough UV-signature check.
    WHY: UV light damages adjacent pyrimidines, so in melanoma most C>T changes
    should have a C or T on the 5' side. A high value confirms the sample is UV-driven.
    CpG share separates the other common C>T source (spontaneous deamination).
    This is a quick check, not a formal signature fit."""
    rows = {}
    for name, c in [("coding", cod), ("noncoding", nc)]:
        total = c.sum()
        ct = c[[x for x in c.index if x[2:5] == "C>T"]]
        dipy = ct[[x for x in ct.index if x[0] in "CT"]].sum()
        cpg = ct[[x for x in ct.index if x[6] == "G"]].sum()
        rows[name] = {"point mutations with context": int(total),
                      "C>T share of all (%)": round(100 * ct.sum() / max(total, 1), 1),
                      "C>T at dipyrimidine, 5' base C or T (%)": round(100 * dipy / max(ct.sum(), 1), 1),
                      "C>T at CpG, 3' base G (%)": round(100 * cpg / max(ct.sum(), 1), 1)}
    t = pd.DataFrame(rows)
    save_table(t, "15_uv_signature_check")
    print(t.to_string())


def table_gene_groups(df):
    """Tables: mutations by GENCODE gene group, and the real content of the 'RNA' class.
    WHY: separates true lncRNA and miRNA mutations from protein-coding genes and
    pseudogenes that the MAF lumps under 'RNA'."""
    order = [g for g in GROUP_ORDER if g in set(df["gene_group"])]
    by_region = pd.crosstab(df["gene_group"], df["region_group"]).reindex(order)
    save_table(by_region, "16_gene_group_by_region")
    rna = df[df["Variant_Classification"] == "RNA"]["gene_group"].value_counts().reindex(order).dropna().astype(int)
    save_table(rna.rename("mutations in the 'RNA' class").to_frame(), "17_rna_class_by_gene_group")
    print("  'RNA' class by GENCODE gene group:")
    print(rna.to_string())
    match = 100 * (df["gene_group"] != "unmatched").mean()
    print(f"  gene symbols matched to GENCODE: {match:.1f}%")
    unmatched = df[df["gene_group"] == "unmatched"]["Hugo_Symbol"].value_counts().head(10)
    if len(unmatched):
        print("  most common unmatched symbols:", ", ".join(f"{k} ({v})" for k, v in unmatched.items()))
    return by_region


def table_true_ncrna(df):
    """Tables: noncoding mutations in true lncRNA and miRNA genes.
    WHY: the first real list of noncoding-RNA mutation candidates. Treat these as
    leads only: no background model has been applied yet."""
    for group, name in [("lncRNA", "18a_lncRNA_genes_mutated"), ("miRNA", "18b_miRNA_genes_mutated")]:
        sub = df[(df["gene_group"] == group) & df["is_noncoding"]]
        t = (sub.groupby("Hugo_Symbol").agg(mutations=("patient", "size"), patients=("patient", "nunique"))
                .sort_values(["patients", "mutations"], ascending=False).head(25))
        save_table(t, name)
        print(f"  noncoding mutations in {group} genes: {len(sub):,} "
              f"({sub['patient'].nunique()} patients, {sub['Hugo_Symbol'].nunique()} genes)")


# ---------------------------------- figures ---------------------------------
def fig_spectrum96(cod, nc):
    """Figure 10. The 96-class mutation spectrum (coding above, noncoding below).
    WHY: the shape of this plot is the mutational process. In melanoma, red C>T
    bars at pyrimidine-pyrimidine contexts should dominate. Similar shapes in
    coding and noncoding DNA justify using one background model for both."""
    fig, axes = plt.subplots(2, 1, figsize=(16, 7.5), sharex=True)
    for ax, counts, title in zip(axes, [cod, nc], ["Coding", "Noncoding"]):
        pct = 100 * counts / max(counts.sum(), 1)
        colors = [SUB_COLORS[c[2:5]] for c in CHANNELS_96]
        ax.bar(range(96), pct.values, color=colors, width=0.8)
        top = max(pct.max(), 1e-9)
        ax.set_ylim(0, top * 1.22)
        for i, sub in enumerate(SPECTRUM):
            ax.add_patch(plt.Rectangle((i * 16 - 0.5, top * 1.12), 16, top * 0.05,
                                       color=SUB_COLORS[sub], clip_on=False))
            ax.text(i * 16 + 7.5, top * 1.19, sub, ha="center", fontsize=9)
        ax.set_ylabel("% of point mutations")
        ax.set_title(f"{title} mutations (n={counts.sum():,})", loc="left", fontsize=10, fontweight="bold", pad=18)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[1].set_xticks(range(96))
    axes[1].set_xticklabels([c[0] + c[2] + c[6] for c in CHANNELS_96], rotation=90, fontsize=6, family="monospace")
    axes[1].set_xlabel("5' base, mutated base, 3' base")
    save_fig(fig, "fig10_spectrum96")


def fig_gene_groups(by_region):
    """Figure 11. Noncoding mutations by GENCODE gene group.
    WHY: shows how many noncoding mutations fall in genuine noncoding-RNA genes
    versus protein-coding genes (introns, UTRs) and pseudogenes."""
    nc_cols = [c for c in by_region.columns if c != "coding"]
    counts = by_region[nc_cols].sum(axis=1)
    counts = counts[counts > 0]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.barh(counts.index, counts.values, color=[GROUP_COLORS.get(g, "#CCCCCC") for g in counts.index])
    ax.invert_yaxis()
    ax.set_xscale("log")
    for y, v in enumerate(counts.values):
        ax.text(v * 1.12, y, f"{v:,}", va="center", fontsize=9)
    ax.set_xlim(counts.min() * 0.6, counts.max() * 6)
    ax.set_xlabel("Noncoding mutations (log scale)")
    ax.set_title("Noncoding mutations by GENCODE gene group", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_fig(fig, "fig11_noncoding_by_gene_group")


# ----------------------------------- main -----------------------------------
def main():
    for d in (TAB, FIG, CTX):
        d.mkdir(parents=True, exist_ok=True)
    src = RES / "all_mutations.csv.gz"
    if not src.exists():
        raise SystemExit(f"Missing {src}. Run 02_merge_and_classify.py first.")

    fasta_path = find_file(FASTA, ["GRCh38.d1.vd1.fa", "GRCh38*.fa", "*.fa", "*.fasta"], "reference FASTA")
    gtf_path = find_file(GTF, ["gencode.v36.annotation.gtf", "gencode*.gtf", "gencode*.gtf.gz"], "GENCODE GTF")
    print(f"FASTA: {fasta_path}\nGTF:   {gtf_path}\n")

    df = pd.read_csv(src, low_memory=False)
    df["Chromosome"] = df["Chromosome"].astype(str)
    print(f"Mutations loaded: {len(df):,}")

    print("\nReading the reference genome (the first run builds an index, about 1-2 minutes)")
    genome = Fasta(str(fasta_path))
    df = add_context(df, genome)
    df = add_spectrum_labels(df)

    print("\nGENCODE gene types")
    genes = load_gene_table(gtf_path)
    df = df.join(genes, on="Hugo_Symbol")
    df["gene_group"] = df["gene_group"].fillna("unmatched")

    print("\nTables")
    table_reference_check(df)
    cod, nc = spectrum96(df)
    table_spectrum96(cod, nc)
    table_uv_check(cod, nc)
    by_region = table_gene_groups(df)
    table_true_ncrna(df)

    print("\nFigures")
    fig_spectrum96(cod, nc)
    fig_gene_groups(by_region)

    out = CTX / "mutations_with_context.csv.gz"
    df.to_csv(out, index=False)
    print(f"\nSaved: {out}\nDone.")


if __name__ == "__main__":
    main()