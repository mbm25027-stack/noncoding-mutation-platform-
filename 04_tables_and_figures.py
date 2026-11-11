"""
04_tables_and_figures.py

Builds every table and figure of the TCGA-SKCM mutation landscape analysis from
the files that scripts 02 and 03 already wrote. It does NOT re-read the 472 MAF
files, so it runs in seconds, and it makes NO PDF:
    tables  -> results/tables/   (CSV, open in Excel)
    figures -> results/figures/  (PNG)

Inputs (created earlier):
    results/all_mutations.csv.gz                    from 02_merge_and_classify.py
    results/coding_explore/coding_mutations.csv.gz  from 03_explore_coding.py
    gdc_sample_sheet*.tsv                           downloaded from the GDC cart

Every function has a WHY note: what the table/figure shows and why the project needs it.
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # write files, never open windows
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ========================= CHANGE THESE IF NEEDED ===========================
BASE = Path(r"D:\noncoding mutation finding surver")   # <<< project folder
GENES = ["BRAF", "NRAS", "TP53"]    # <<< genes for position tables and lollipop plots
HIGHLIGHT = {"BRAF", "NRAS"}        # <<< genes drawn in the highlight colour in figures 6 and 7
TOP_N = 25                          # <<< how many rows/bars in "top" tables and figures
MIN_PATIENTS_NONCODING = 3          # <<< recurrent noncoding positions need at least this many patients
META = {                            # <<< provenance, copied into the dataset summary table
    "Project": "TCGA-SKCM",
    "GDC data release": "46.0",
    "Experimental strategy": "WXS (whole-exome sequencing)",
    "Data type": "Masked Somatic Mutation (open access, MAF)",
    "Download date": "2026-11-07",
}
# ============================================================================

RES = BASE / "results"
TAB = RES / "tables"
FIG = RES / "figures"

GROUP_ORDER = ["coding", "splice", "intron", "ncRNA", "UTR", "flank", "intergenic", "other"]
GROUP_COLORS = {
    "coding": "#1F3864", "splice": "#2A9D8F", "intron": "#E76F51", "ncRNA": "#6A4C93",
    "UTR": "#F4A261", "flank": "#8D99AE", "intergenic": "#999999", "other": "#CCCCCC",
}
CLASS_COLORS = {
    "Missense_Mutation": "tab:green", "Nonsense_Mutation": "tab:red",
    "Frame_Shift_Del": "tab:purple", "Frame_Shift_Ins": "tab:purple",
    "In_Frame_Del": "tab:orange", "In_Frame_Ins": "tab:orange",
    "Silent": "lightgray", "Nonstop_Mutation": "tab:brown",
    "Translation_Start_Site": "tab:pink",
}
SAMPLE_TYPE = {"01": "Primary tumour (01)", "06": "Metastatic (06)", "02": "Recurrent (02)",
               "03": "Blood-derived (03)", "05": "Additional primary (05)"}
SPECTRUM = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BLUE, ORANGE = "#1F3864", "#E76F51"


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


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def pct(n, total):
    return round(100 * n / total, 2)


def substitution_spectrum(df):
    """Count point mutations per substitution class, folded onto the pyrimidine
    (C or T) strand so that G>A is counted as C>T and A>G as T>C."""
    snv = df[(df["Variant_Type"] == "SNP")
             & df["Reference_Allele"].isin(list("ACGT"))
             & df["Tumor_Seq_Allele2"].isin(list("ACGT"))]
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    flip = snv["Reference_Allele"].isin(["A", "G"])
    ref = snv["Reference_Allele"].where(~flip, snv["Reference_Allele"].map(comp))
    alt = snv["Tumor_Seq_Allele2"].where(~flip, snv["Tumor_Seq_Allele2"].map(comp))
    return (ref + ">" + alt).value_counts().reindex(SPECTRUM).fillna(0).astype(int)


# --------------------------------- loading ----------------------------------
def load_data():
    """Load the two tables written by scripts 02 and 03 and check they agree.

    WHY: starting from the saved tables means every number in the tables and
    figures comes from one cleaned dataset (one tumour sample per patient,
    de-duplicated), so figures can never disagree with each other."""
    all_path = RES / "all_mutations.csv.gz"
    cod_path = RES / "coding_explore" / "coding_mutations.csv.gz"
    for p, script in [(all_path, "02_merge_and_classify.py"), (cod_path, "03_explore_coding.py")]:
        if not p.exists():
            raise SystemExit(f"Missing {p}. Run {script} first.")
    allm = pd.read_csv(all_path, low_memory=False)
    coding = pd.read_csv(cod_path, low_memory=False)

    need_all = {"patient", "sample_code", "Variant_Classification", "Variant_Type", "region_group",
                "is_noncoding", "Hugo_Symbol", "Chromosome", "Start_Position",
                "Reference_Allele", "Tumor_Seq_Allele2", "NCBI_Build"}
    need_cod = {"patient", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short", "Protein_position"}
    for name, df, need in [("all_mutations", allm, need_all), ("coding_mutations", coding, need_cod)]:
        missing = need - set(df.columns)
        if missing:
            raise SystemExit(f"{name} is missing columns {sorted(missing)}. Re-run scripts 02 and 03.")

    n_cod_all = int((allm["region_group"] == "coding").sum())
    if n_cod_all != len(coding):
        print(f"WARNING: {n_cod_all:,} coding rows in all_mutations but {len(coding):,} in "
              "coding_mutations. Re-run 02 and 03 together so both use the same data.")
    return allm, coding


# ---------------------------------- tables ----------------------------------
def table_dataset_summary(allm, per_patient):
    """Table: dataset summary.
    WHY: one place that records exactly what data the analysis used (release,
    filters, counts). It becomes the first lines of your methods section."""
    n_files = len(list((BASE / "maf_files").rglob("*.maf.gz")))
    n_nc = int(allm["is_noncoding"].sum())
    types = allm.drop_duplicates("patient")["sample_code"].astype(str).str.zfill(2)
    rows = list(META.items()) + [
        ("MAF files analysed", f"{n_files}"),
        ("Patients (one tumour sample each)", f"{allm['patient'].nunique()}"),
        ("Mutations after de-duplication", f"{len(allm):,}"),
        ("Noncoding mutations", f"{n_nc:,} ({pct(n_nc, len(allm))}%)"),
        ("Median mutations per patient", f"{per_patient.median():.0f}"),
        ("Genome build", ", ".join(sorted(allm['NCBI_Build'].dropna().astype(str).unique()))),
    ]
    for code, n in types.value_counts().items():
        rows.append((f"Samples: {SAMPLE_TYPE.get(code, code)}", f"{n}"))
    save_table(pd.DataFrame(rows, columns=["Item", "Value"]), "01_dataset_summary", index=False)


def table_variant_types(allm):
    """Table: SNP / insertion / deletion counts.
    WHY: shows the data are almost all point mutations, so SNV-based methods
    (substitution spectrum, background model) cover nearly everything."""
    t = allm["Variant_Type"].value_counts().rename("count").to_frame()
    t["percent"] = (100 * t["count"] / t["count"].sum()).round(2)
    save_table(t, "02_variant_types")


def table_variant_classes(allm):
    """Table: mutations per Variant_Classification and region group.
    WHY: this is the step that separates coding from noncoding mutations, the
    basis of the whole noncoding project."""
    t = allm["Variant_Classification"].value_counts().rename("mutations").to_frame()
    t["percent"] = (100 * t["mutations"] / t["mutations"].sum()).round(2)
    grp = allm.drop_duplicates("Variant_Classification").set_index("Variant_Classification")["region_group"]
    t["region_group"] = t.index.map(grp)
    save_table(t, "03_variant_classification")
    return t


def table_region_groups(allm):
    """Table: mutations per region group (coding, splice, intron, ncRNA, UTR, flank, intergenic).
    WHY: quantifies how little noncoding signal an exome holds, which is the
    reason the final analysis needs whole-genome data."""
    g = allm["region_group"].value_counts()
    t = g.rename("mutations").to_frame()
    t["percent"] = (100 * t["mutations"] / t["mutations"].sum()).round(2)
    save_table(t, "04_region_groups")
    return g


def table_burden(per_patient):
    """Table: mutations per patient (summary and the full list).
    WHY: mutation burden varies enormously between patients. A few hypermutated
    tumours can create fake 'recurrent' positions, so the background model must
    use per-patient rates."""
    q = per_patient.quantile([0, .05, .25, .5, .75, .95, 1])
    s = pd.DataFrame({"mutations_per_patient": q.round(0).astype(int)})
    s.index = ["minimum", "5th percentile", "25th percentile", "median",
               "75th percentile", "95th percentile", "maximum"]
    save_table(s, "05_mutation_burden_summary")
    save_table(per_patient.rename("mutations").to_frame(), "05b_mutations_per_patient")


def table_spectrum(coding_counts, nc_counts):
    """Table: substitution spectrum for coding and noncoding point mutations.
    WHY: the pattern of base changes is a fingerprint of the mutational process.
    In melanoma, C>T points to UV damage. The background model must reproduce
    this pattern, otherwise it will call UV-prone sites 'hotspots'."""
    t = pd.DataFrame({
        "coding_count": coding_counts, "coding_percent": (100 * coding_counts / coding_counts.sum()).round(1),
        "noncoding_count": nc_counts, "noncoding_percent": (100 * nc_counts / nc_counts.sum()).round(1),
    })
    save_table(t, "06_substitution_spectrum")


def table_recurrent_protein(coding):
    """Table: most recurrent protein-level changes (non-silent).
    WHY: positive control. If BRAF V600E, NRAS Q61 and other known melanoma
    hotspots come out on top, the data and the pipeline are behaving correctly."""
    n_pat = coding["patient"].nunique()
    prot = coding[(coding["Variant_Classification"] != "Silent") & coding["HGVSp_Short"].notna()]
    r = prot.groupby(["Hugo_Symbol", "HGVSp_Short"])["patient"].nunique().sort_values(ascending=False).head(TOP_N)
    t = r.rename("patients").to_frame()
    t["percent"] = (100 * t["patients"] / n_pat).round(1)
    save_table(t, "07_recurrent_protein_changes")
    return t


def table_top_genes(coding):
    """Table: genes with non-silent coding mutations in the most patients.
    WHY: shows the size bias. Huge genes (TTN, MUC16) top the list because they
    offer more places to mutate, not because they drive cancer. This is the
    same trap that makes raw recurrence misleading in noncoding regions."""
    n_pat = coding["patient"].nunique()
    ns = coding[coding["Variant_Classification"] != "Silent"]
    g = ns.groupby("Hugo_Symbol")["patient"].nunique().sort_values(ascending=False).head(TOP_N)
    t = g.rename("patients").to_frame()
    t["percent"] = (100 * t["patients"] / n_pat).round(1)
    save_table(t, "08_top_mutated_genes")
    return t


def gene_positions(coding, gene):
    """Return coding mutations of one gene with a numeric protein position."""
    sub = coding[coding["Hugo_Symbol"] == gene].copy()
    sub["pos"] = pd.to_numeric(sub["Protein_position"].astype(str).str.extract(r"(\d+)")[0], errors="coerce")
    return sub.dropna(subset=["pos"])


def table_driver_positions(coding):
    """Tables: mutation position and class for chosen genes (BRAF, NRAS, TP53).
    WHY: answers 'which kind of mutation occurs at which position'. Known
    drivers show a few sharp peaks (hotspots), which is what a hotspot caller
    should detect."""
    out = {}
    for gene in GENES:
        sub = gene_positions(coding, gene)
        if sub.empty:
            print(f"  note: no mutations with a protein position for {gene}")
            continue
        g = (sub.groupby(["pos", "Variant_Classification"])["patient"].nunique()
                .reset_index(name="patients").sort_values("patients", ascending=False))
        g["pos"] = g["pos"].astype(int)
        save_table(g.head(15), f"09_positions_{gene}", index=False)
        out[gene] = (sub, g)
    return out


def table_ncrna_genes(allm):
    """Table: genes with most mutations in the 'RNA' class.
    WHY: first look at noncoding-RNA mutations. It also exposes a labelling
    problem: protein-coding genes (SSPOP, DCHS2) and pseudogenes appear here,
    so true lncRNA/miRNA must later be defined from GENCODE gene types."""
    r = allm[allm["Variant_Classification"] == "RNA"]["Hugo_Symbol"].value_counts().head(15)
    t = r.rename("mutations").to_frame()
    save_table(t, "10_ncrna_class_genes")
    return t


def table_noncoding_recurrent(allm):
    """Table: noncoding positions mutated in several patients.
    WHY: a preview of hotspot calling, and a warning. Most recurrent positions
    are C>T at UV-prone sites or in huge/repetitive genes (TTN, pseudogenes),
    so recurrence alone is not evidence of a hotspot."""
    nc = allm[allm["is_noncoding"]]
    keys = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
            "Variant_Classification", "Hugo_Symbol"]
    r = nc.groupby(keys)["patient"].nunique().sort_values(ascending=False)
    r = r[r >= MIN_PATIENTS_NONCODING].rename("patients").reset_index()
    save_table(r, "11_noncoding_recurrent_positions", index=False)
    return r


def table_missing_cases(allm):
    """Table: cases listed in the GDC sample sheet but absent from the analysis.
    WHY: the sample sheet lists 470 cases but 467 patients were analysed. Finding
    the missing 3 (empty or unreadable files?) protects the integrity of every
    later step."""
    sheets = sorted(BASE.glob("gdc_sample_sheet*.tsv"))
    if not sheets:
        print("  note: no gdc_sample_sheet*.tsv in the project folder, reconciliation skipped")
        return
    ss = pd.read_csv(sheets[0], sep="\t", dtype=str)
    case_col = next((c for c in ss.columns if c.strip().lower() == "case id"), None)
    if case_col is None:
        print(f"  note: no 'Case ID' column in the sample sheet (columns: {list(ss.columns)})")
        return
    pat = re.compile(r"TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
    ss["_cases"] = ss[case_col].fillna("").apply(pat.findall)
    expected = {c for cs in ss["_cases"] for c in cs}
    observed = set(allm["patient"])
    missing = sorted(expected - observed)
    print(f"  reconciliation: {len(expected)} cases in sample sheet, {len(observed)} patients analysed, "
          f"{len(missing)} missing")
    save_table(pd.DataFrame({"case": missing}), "12_missing_cases", index=False)
    rows = ss[ss["_cases"].apply(lambda cs: any(c in missing for c in cs))].drop(columns="_cases")
    save_table(rows, "12b_missing_cases_sample_sheet_rows", index=False)


# ---------------------------------- figures ---------------------------------
def fig_region_groups(group_counts):
    """Figure 1. Mutations by region group (log scale).
    WHY: one glance shows that coding mutations outnumber every noncoding
    class by orders of magnitude in exome data. Log scale keeps the small
    classes visible."""
    order = [g for g in GROUP_ORDER if g in group_counts.index]
    vals = group_counts.reindex(order)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.barh(order, vals.values, color=[GROUP_COLORS[g] for g in order])
    ax.set_xscale("log")
    ax.set_xlim(1, vals.max() * 8)
    ax.invert_yaxis()
    for y, v in enumerate(vals.values):
        ax.text(v * 1.15, y, f"{v:,}", va="center", fontsize=9)
    ax.set_xlabel("Mutations (log scale)")
    ax.set_title("Mutations by genomic region group", fontweight="bold")
    clean_axes(ax)
    save_fig(fig, "fig1_region_groups")


def fig_coding_vs_noncoding(group_counts):
    """Figure 2. Coding versus noncoding share.
    WHY: the same message as figure 1 on a linear scale, so the real size of
    the noncoding fraction is not visually exaggerated by the log axis."""
    n_nc = int(group_counts.drop(labels=["coding"], errors="ignore").sum())
    n_cod = int(group_counts.get("coding", 0))
    total = n_cod + n_nc
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.bar(["Coding", "Noncoding"], [n_cod, n_nc], color=[BLUE, "#2A9D8F"])
    for i, v in enumerate([n_cod, n_nc]):
        ax.text(i, v, f"{100 * v / total:.1f}%\n({v:,})", ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Mutations")
    ax.set_ylim(0, n_cod * 1.18)
    ax.set_title("Coding versus noncoding share", fontweight="bold")
    clean_axes(ax)
    save_fig(fig, "fig2_coding_vs_noncoding")


def fig_variant_classes(classes, top=12):
    """Figure 3. Top variant classifications, coloured by region group.
    WHY: shows which specific kinds of mutation dominate (missense, silent,
    nonsense) and which noncoding kinds exist (intron, splice, UTR, flank)."""
    t = classes.head(top)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.barh(t.index, t["mutations"], color=[GROUP_COLORS.get(g, "#CCCCCC") for g in t["region_group"]])
    ax.set_xscale("log")
    ax.set_xlim(t["mutations"].min() * 0.5, t["mutations"].max() * 8)
    ax.invert_yaxis()
    for y, (n, p) in enumerate(zip(t["mutations"], t["percent"])):
        ax.text(n * 1.15, y, f"{n:,} ({p:.2f}%)", va="center", fontsize=8)
    handles = [plt.Rectangle((0, 0), 1, 1, color=GROUP_COLORS[g]) for g in GROUP_ORDER if g in set(t["region_group"])]
    ax.legend(handles, [g for g in GROUP_ORDER if g in set(t["region_group"])], title="Region group",
              frameon=False, loc="lower right", fontsize=8)
    ax.set_xlabel("Mutations (log scale)")
    ax.set_title(f"Top {top} variant classifications", fontweight="bold")
    clean_axes(ax)
    save_fig(fig, "fig3_variant_classes")


def fig_burden(per_patient):
    """Figure 4. Mutations per patient.
    WHY: reveals hypermutated tumours and the wide spread of burden. Uses
    log-spaced bins (equal-width bins on a log axis would distort the shape)."""
    v = per_patient.values
    bins = np.logspace(np.log10(v.min()), np.log10(v.max()), 30)
    med = float(np.median(v))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(v, bins=bins, color=BLUE)
    ax.set_xscale("log")
    ax.axvline(med, color=ORANGE, lw=2)
    ax.text(med * 1.08, ax.get_ylim()[1] * 0.92, f"median {med:.0f}", color=ORANGE)
    ax.set_xlabel("Mutations per patient (log scale)")
    ax.set_ylabel("Patients")
    ax.set_title("Mutational burden per patient", fontweight="bold")
    clean_axes(ax)
    save_fig(fig, "fig4_mutation_burden")


def fig_spectrum(coding_counts, nc_counts):
    """Figure 5. Substitution spectrum, coding versus noncoding.
    WHY: C>T dominance is the UV signature of melanoma. Showing noncoding
    mutations next to coding ones checks that the same process acts in both,
    which is the assumption behind using one background model."""
    x = np.arange(len(SPECTRUM))
    w = 0.38
    cod = 100 * coding_counts / coding_counts.sum()
    nc = 100 * nc_counts / nc_counts.sum() if nc_counts.sum() else nc_counts * 0
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(x - w / 2, cod.values, w, color=BLUE, label=f"Coding (n={coding_counts.sum():,})")
    ax.bar(x + w / 2, nc.values, w, color=ORANGE, label=f"Noncoding (n={nc_counts.sum():,})")
    for xi, v in zip(x - w / 2, cod.values):
        ax.text(xi, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    for xi, v in zip(x + w / 2, nc.values):
        ax.text(xi, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(SPECTRUM)
    ax.set_ylabel("% of point mutations")
    ax.set_title("Substitution spectrum", fontweight="bold")
    ax.legend(frameon=False)
    clean_axes(ax)
    save_fig(fig, "fig5_substitution_spectrum")


def fig_recurrent_protein(t, n_patients):
    """Figure 6. Most recurrent protein-level changes.
    WHY: visual positive control. Known melanoma hotspots should dominate."""
    labels = [f"{g} {p}" for g, p in t.index]
    colors = [ORANGE if g in HIGHLIGHT else BLUE for g, _ in t.index]
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(labels, t["percent"], color=colors)
    ax.invert_yaxis()
    for y, (n, p) in enumerate(zip(t["patients"], t["percent"])):
        ax.text(p + 0.4, y, f"{n} pts ({p:.1f}%)", va="center", fontsize=8)
    ax.set_xlim(0, t["percent"].max() * 1.25)
    ax.set_xlabel(f"% of {n_patients} patients")
    ax.set_title("Most recurrent protein-level changes", fontweight="bold")
    clean_axes(ax)
    save_fig(fig, "fig6_recurrent_protein_changes")


def fig_top_genes(t):
    """Figure 7. Genes mutated in the most patients (size-biased list).
    WHY: demonstrates why a mutation count per gene is not a driver test and
    why a background model is required."""
    colors = [ORANGE if g in HIGHLIGHT else "#8D99AE" for g in t.index]
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.barh(t.index, t["percent"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("% of patients with a non-silent coding mutation")
    ax.set_title("Genes mutated in the most patients (size-biased list)", fontweight="bold")
    clean_axes(ax)
    save_fig(fig, "fig7_top_mutated_genes")


def fig_lollipop(gene, sub, n_patients):
    """Figure 8. Mutation position along one protein (lollipop plot).
    WHY: shows 'what kind of mutation at which position'. Drivers such as BRAF
    show one tall peak (a hotspot); tumour suppressors such as TP53 show many
    different mutations spread along the protein."""
    length = pd.to_numeric(sub["Protein_position"].astype(str).str.extract(r"/(\d+)")[0], errors="coerce").max()
    g = sub.groupby(["pos", "Variant_Classification"])["patient"].nunique().reset_index(name="patients")
    n_gene = sub[sub["Variant_Classification"] != "Silent"]["patient"].nunique()
    fig, ax = plt.subplots(figsize=(11, 3.8))
    for c, part in g.groupby("Variant_Classification"):
        color = CLASS_COLORS.get(c, "black")
        ax.vlines(part["pos"], 0, part["patients"], color=color, lw=1)
        ax.scatter(part["pos"], part["patients"], s=28, color=color, label=c, zorder=3)
    if pd.notna(length):
        ax.set_xlim(0, length)
    ax.set_xlabel("Position along the protein (amino acid number)")
    ax.set_ylabel("Patients")
    ax.set_title(f"{gene}: coding mutations, {n_gene} of {n_patients} patients ({100 * n_gene / n_patients:.1f}%) "
                 "with a non-silent mutation", fontweight="bold", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    clean_axes(ax)
    save_fig(fig, f"fig8_lollipop_{gene}")


def fig_ncrna(t):
    """Figure 9. Genes with most mutations in the 'RNA' class.
    WHY: first look at noncoding-RNA mutations, with the caution that many
    entries are pseudogenes or protein-coding genes and prone to mapping errors."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(t.index, t["mutations"], color=GROUP_COLORS["ncRNA"])
    ax.invert_yaxis()
    ax.set_xlabel("Mutations in the 'RNA' class")
    ax.set_title("Most frequently mutated 'RNA'-class genes", fontweight="bold")
    clean_axes(ax)
    save_fig(fig, "fig9_ncrna_class_genes")


# ----------------------------------- main -----------------------------------
def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    allm, coding = load_data()
    n_patients = allm["patient"].nunique()
    per_patient = allm.groupby("patient").size().sort_values()
    print(f"Patients: {n_patients}   Mutations: {len(allm):,}\n")

    print("Tables")
    table_dataset_summary(allm, per_patient)
    table_variant_types(allm)
    classes = table_variant_classes(allm)
    groups = table_region_groups(allm)
    table_burden(per_patient)
    cod_spec = substitution_spectrum(allm[allm["region_group"] == "coding"])
    nc_spec = substitution_spectrum(allm[allm["is_noncoding"]])
    table_spectrum(cod_spec, nc_spec)
    recur = table_recurrent_protein(coding)
    genes_t = table_top_genes(coding)
    drivers = table_driver_positions(coding)
    ncrna = table_ncrna_genes(allm)
    table_noncoding_recurrent(allm)
    table_missing_cases(allm)

    print("\nFigures")
    fig_region_groups(groups)
    fig_coding_vs_noncoding(groups)
    fig_variant_classes(classes)
    fig_burden(per_patient)
    fig_spectrum(cod_spec, nc_spec)
    fig_recurrent_protein(recur, coding["patient"].nunique())
    fig_top_genes(genes_t)
    for gene, (sub, _) in drivers.items():
        fig_lollipop(gene, sub, coding["patient"].nunique())
    fig_ncrna(ncrna)

    print(f"\nDone. Tables: {TAB}\n      Figures: {FIG}")


if __name__ == "__main__":
    main()