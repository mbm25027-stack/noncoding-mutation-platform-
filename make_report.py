"""
make_report.py  -  build the melanoma mutation-landscape PDF report from raw GDC MAF files.

USAGE (PowerShell, from your project folder):
    pip install pandas matplotlib reportlab
    python make_report.py
or with your own paths:
    python make_report.py --maf-dir "D:\\noncoding mutation finding surver\\maf_files" --out-dir "D:\\noncoding mutation finding surver\\report"

OUTPUT (in --out-dir):
    Mutation_Landscape_Report.pdf     the report
    figures/*.png                     every chart as a separate image
    tables/*.csv                      every table as a CSV (for Excel / your own plots)

What it does:
    1. reads every .maf / .maf.gz under --maf-dir
    2. keeps ONE tumour sample per patient (prefers primary '01', then metastatic '06')
    3. removes duplicate mutation rows
    4. classifies mutations into coding / splice / intron / ncRNA / UTR / flank / intergenic
    5. computes statistics, draws charts, writes the PDF
Nothing is hard-coded: every number in the report is computed from your files.
"""
import argparse
import glob
import html
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (Image, KeepTogether, PageBreak, Paragraph,
                                SimpleDocTemplate, Spacer, Table, TableStyle)

# ----------------------------------------------------------------------------
# Settings you may want to change
# ----------------------------------------------------------------------------
DEFAULT_MAF_DIR = r"D:\noncoding mutation finding surver\maf_files"
DEFAULT_OUT_DIR = r"D:\noncoding mutation finding surver\report"
REPORT_TITLE = "Somatic Mutation Landscape of Melanoma"
REPORT_SUBTITLE = "Coding and noncoding mutation survey from GDC open-access MAF files"
DRIVER_GENES = ["BRAF", "NRAS", "TP53"]      # genes to examine position by position
MIN_PATIENTS_NONCODING = 3                    # recurrence cut-off for noncoding positions

# Variant_Classification -> region group
REGION = {
    "Missense_Mutation": "coding", "Silent": "coding", "Nonsense_Mutation": "coding",
    "Frame_Shift_Del": "coding", "Frame_Shift_Ins": "coding", "In_Frame_Del": "coding",
    "In_Frame_Ins": "coding", "Translation_Start_Site": "coding", "Nonstop_Mutation": "coding",
    "Splice_Site": "splice", "Splice_Region": "splice",
    "Intron": "intron", "RNA": "ncRNA",
    "5'UTR": "UTR", "3'UTR": "UTR", "5'Flank": "flank", "3'Flank": "flank", "IGR": "intergenic",
}
SAMPLE_RANK = {"01": 0, "06": 1}              # lower = preferred when a patient has several samples
SAMPLE_NAME = {"01": "Primary tumour (01)", "06": "Metastatic (06)", "02": "Recurrent (02)",
               "03": "Blood-derived (03)", "05": "Additional new primary (05)"}

WANTED = ["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
          "Variant_Classification", "Variant_Type", "NCBI_Build", "Tumor_Sample_Barcode", "HGVSp_Short"]

NAVY, TEAL, ORANGE, GREY, PURPLE, SAND = "#1F3A5F", "#2A9D8F", "#E76F51", "#8D99AE", "#6A4C93", "#F4A261"
GROUP_COLOURS = {"coding": NAVY, "splice": TEAL, "intron": ORANGE, "ncRNA": PURPLE,
                 "UTR": SAND, "flank": GREY, "intergenic": "#999999", "other": "#BBBBBB"}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


# ----------------------------------------------------------------------------
# 1. Load and clean
# ----------------------------------------------------------------------------
def find_maf_files(maf_dir):
    files = []
    for pattern in ("*.maf", "*.maf.gz", "*.maf.txt", "*.maf.txt.gz"):
        files += glob.glob(os.path.join(maf_dir, "**", pattern), recursive=True)
    return sorted(set(files))


def load_mafs(files):
    frames = []
    for i, f in enumerate(files, 1):
        # comment="#" skips the '#version' header lines GDC puts at the top of each MAF
        d = pd.read_csv(f, sep="\t", comment="#", dtype=str, low_memory=False,
                        usecols=lambda c: c in WANTED)
        frames.append(d)
        if i % 100 == 0:
            print(f"  read {i}/{len(files)} files")
    return pd.concat(frames, ignore_index=True)


def clean(df):
    info = {"rows_read": len(df)}
    df = df.dropna(subset=["Tumor_Sample_Barcode", "Variant_Classification"]).copy()
    if "HGVSp_Short" not in df.columns:
        df["HGVSp_Short"] = pd.NA
    df["patient"] = df["Tumor_Sample_Barcode"].str[:12]
    df["sample_code"] = df["Tumor_Sample_Barcode"].str[13:15]

    # keep one tumour sample per patient
    samples = df[["patient", "Tumor_Sample_Barcode", "sample_code"]].drop_duplicates()
    info["multi_sample_patients"] = int((samples.groupby("patient").size() > 1).sum())
    samples["rank"] = samples["sample_code"].map(SAMPLE_RANK).fillna(2)
    samples = samples.sort_values(["patient", "rank", "Tumor_Sample_Barcode"]).drop_duplicates("patient")
    df = df[df["Tumor_Sample_Barcode"].isin(set(samples["Tumor_Sample_Barcode"]))]
    info["sample_types"] = samples["sample_code"].value_counts()

    # remove duplicate mutation rows within a patient
    df = df.drop_duplicates(subset=["patient", "Chromosome", "Start_Position",
                                    "Reference_Allele", "Tumor_Seq_Allele2"]).copy()
    df["region_group"] = df["Variant_Classification"].map(REGION).fillna("other")
    info["n_patients"] = df["patient"].nunique()
    info["n_mutations"] = len(df)
    return df, info


# ----------------------------------------------------------------------------
# 2. Analysis (every function returns a DataFrame or a number)
# ----------------------------------------------------------------------------
def pyrimidine_spectrum(coding):
    """Substitution spectrum of coding SNPs, folded onto the pyrimidine (C/T) strand."""
    snp = coding[(coding["Variant_Type"] == "SNP") &
                 coding["Reference_Allele"].isin(COMPLEMENT) &
                 coding["Tumor_Seq_Allele2"].isin(COMPLEMENT)]
    ref, alt = snp["Reference_Allele"], snp["Tumor_Seq_Allele2"]
    flip = ref.isin(["G", "A"])
    ref = ref.where(~flip, ref.map(COMPLEMENT))
    alt = alt.where(~flip, alt.map(COMPLEMENT))
    s = (ref + ">" + alt).value_counts()
    order = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    s = s.reindex(order).fillna(0).astype(int)
    out = pd.DataFrame({"count": s})
    out["percent"] = (100 * out["count"] / max(out["count"].sum(), 1)).round(1)
    return out


def top_recurrent_protein(nonsilent, n_patients, top=15):
    d = nonsilent.dropna(subset=["HGVSp_Short"])
    d = d[d["HGVSp_Short"] != ""]
    t = d.groupby(["Hugo_Symbol", "HGVSp_Short"])["patient"].nunique().sort_values(ascending=False).head(top)
    t = t.reset_index().rename(columns={"patient": "patients"})
    t["percent"] = (100 * t["patients"] / n_patients).round(1)
    return t


def top_genes(nonsilent, n_patients, top=25):
    t = nonsilent.groupby("Hugo_Symbol")["patient"].nunique().sort_values(ascending=False).head(top)
    t = t.reset_index().rename(columns={"patient": "patients"})
    t["percent"] = (100 * t["patients"] / n_patients).round(1)
    return t


def gene_positions(coding, gene, top=8):
    d = coding[(coding["Hugo_Symbol"] == gene)].copy()
    d["pos"] = d["HGVSp_Short"].str.extract(r"p\.[A-Za-z\*]+(\d+)", expand=False)
    d = d.dropna(subset=["pos"])
    t = (d.groupby(["pos", "Variant_Classification"])["patient"].nunique()
           .sort_values(ascending=False).head(top).reset_index())
    return t.rename(columns={"patient": "patients"})


def noncoding_recurrent(noncoding, min_patients, top=20):
    keys = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
            "Variant_Classification", "Hugo_Symbol"]
    t = noncoding.groupby(keys)["patient"].nunique()
    t = t[t >= min_patients].sort_values(ascending=False).head(top).reset_index()
    t["position"] = (t["Chromosome"] + ":" + t["Start_Position"] + " " +
                     t["Reference_Allele"] + ">" + t["Tumor_Seq_Allele2"])
    return t.rename(columns={"patient": "patients"})[["position", "Variant_Classification",
                                                      "Hugo_Symbol", "patients"]]


# ----------------------------------------------------------------------------
# 3. Figures
# ----------------------------------------------------------------------------
def save_fig(fig, path):
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path


def make_figures(df, res, fig_dir):
    F = {}
    # region groups
    rg = res["region_counts"].sort_values()
    fig, ax = plt.subplots(figsize=(7, 3.4))
    ax.barh(rg.index, rg.values, color=[GROUP_COLOURS.get(g, GREY) for g in rg.index])
    ax.set_xscale("log"); ax.set_xlabel("Mutations (log scale)")
    for i, v in enumerate(rg.values):
        ax.text(v * 1.12, i, f"{v:,}", va="center", fontsize=9)
    ax.set_xlim(1, rg.max() * 12)
    ax.set_title("Mutations by genomic region group", loc="left", fontweight="bold")
    F["region"] = save_fig(fig, os.path.join(fig_dir, "fig1_region_groups.png"))

    # coding vs noncoding
    n_cod, n_non = res["n_coding"], res["n_noncoding"]
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.bar([f"Coding\n({n_cod:,})", f"Noncoding\n({n_non:,})"], [n_cod, n_non], color=[NAVY, TEAL])
    tot = n_cod + n_non
    for i, v in enumerate([n_cod, n_non]):
        ax.text(i, v, f"{100 * v / tot:.1f}%", ha="center", va="bottom")
    ax.set_ylabel("Mutations")
    ax.set_title("Coding versus noncoding share", loc="left", fontweight="bold")
    F["share"] = save_fig(fig, os.path.join(fig_dir, "fig2_coding_vs_noncoding.png"))

    # classifications
    vc = res["class_table"].head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.barh(vc["classification"], vc["mutations"], color=[GROUP_COLOURS.get(g, GREY) for g in vc["region_group"]])
    ax.set_xscale("log"); ax.set_xlabel("Mutations (log scale)")
    for i, (m, p) in enumerate(zip(vc["mutations"], vc["percent"])):
        ax.text(m * 1.1, i, f"{m:,}  ({p:.2f}%)", va="center", fontsize=8)
    ax.set_xlim(vc["mutations"].min() / 2, vc["mutations"].max() * 12)
    ax.set_title("Top variant classifications (colour = region group)", loc="left", fontweight="bold")
    F["classes"] = save_fig(fig, os.path.join(fig_dir, "fig3_variant_classes.png"))

    # substitution spectrum
    sp = res["spectrum"]
    fig, ax = plt.subplots(figsize=(7, 3.3))
    top_sub = sp["percent"].idxmax()
    ax.bar(sp.index, sp["percent"], color=[ORANGE if s == top_sub else GREY for s in sp.index])
    for i, v in enumerate(sp["percent"]):
        ax.text(i, v + 1.2, f"{v}%", ha="center", fontsize=9)
    ax.set_ylabel("% of coding point mutations"); ax.set_ylim(0, 100)
    ax.set_title(f"Substitution spectrum ({top_sub} dominates)", loc="left", fontweight="bold")
    F["spectrum"] = save_fig(fig, os.path.join(fig_dir, "fig4_substitution_spectrum.png"))

    # recurrent protein changes
    rp = res["recurrent"].iloc[::-1]
    labels = rp["Hugo_Symbol"] + " " + rp["HGVSp_Short"]
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.barh(labels, rp["percent"], color=[ORANGE if g in ("BRAF", "NRAS") else NAVY for g in rp["Hugo_Symbol"]])
    for i, (n, p) in enumerate(zip(rp["patients"], rp["percent"])):
        ax.text(p + 0.4, i, f"{n} pts ({p}%)", va="center", fontsize=8)
    ax.set_xlabel(f"% of {res['n_patients']} patients"); ax.set_xlim(0, max(rp["percent"].max() * 1.3, 5))
    ax.set_title("Most recurrent protein-level changes", loc="left", fontweight="bold")
    F["recurrent"] = save_fig(fig, os.path.join(fig_dir, "fig5_recurrent_protein_changes.png"))

    # top genes
    tg = res["top_genes"].iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.barh(tg["Hugo_Symbol"], tg["percent"], color=[ORANGE if g in DRIVER_GENES else GREY for g in tg["Hugo_Symbol"]])
    ax.set_xlabel("% of patients with a non-silent coding mutation")
    ax.set_title("Genes mutated in the most patients (size-biased list)", loc="left", fontweight="bold")
    F["genes"] = save_fig(fig, os.path.join(fig_dir, "fig6_top_genes.png"))

    # ncRNA genes
    nr = res["ncrna_genes"].iloc[::-1]
    if len(nr):
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.barh(nr["Hugo_Symbol"], nr["mutations"], color=TEAL)
        ax.set_xlabel("Mutations in the 'RNA' class")
        ax.set_title("Most frequently mutated ncRNA-class genes", loc="left", fontweight="bold")
        F["ncrna"] = save_fig(fig, os.path.join(fig_dir, "fig7_ncrna_genes.png"))

    # mutations per patient
    pp = res["per_patient"]
    fig, ax = plt.subplots(figsize=(7, 3.2))
    bins = pd.Series(pp).clip(lower=1)
    ax.hist(bins, bins=40, color=NAVY)
    ax.set_xscale("log"); ax.set_xlabel("Mutations per patient (log scale)"); ax.set_ylabel("Patients")
    ax.axvline(pp.median(), color=ORANGE, lw=1.5)
    ax.text(pp.median() * 1.05, ax.get_ylim()[1] * 0.9, f"median {int(pp.median()):,}", color=ORANGE, fontsize=9)
    ax.set_title("Mutational burden per patient", loc="left", fontweight="bold")
    F["burden"] = save_fig(fig, os.path.join(fig_dir, "fig8_mutations_per_patient.png"))
    return F


# ----------------------------------------------------------------------------
# 4. PDF
# ----------------------------------------------------------------------------
def esc(x):
    return html.escape(str(x), quote=False)


def build_pdf(pdf_path, res, info, F):
    ss = getSampleStyleSheet()
    H1 = ParagraphStyle("H1", parent=ss["Heading1"], textColor=colors.HexColor(NAVY), fontSize=15,
                        spaceBefore=10, spaceAfter=6, keepWithNext=1)
    H2 = ParagraphStyle("H2", parent=ss["Heading2"], textColor=colors.HexColor(NAVY), fontSize=11.5,
                        spaceBefore=8, spaceAfter=4, keepWithNext=1)
    B = ParagraphStyle("B", parent=ss["BodyText"], fontSize=9.8, leading=14, alignment=TA_JUSTIFY, spaceAfter=5)
    BL = ParagraphStyle("BL", parent=B, leftIndent=14, bulletIndent=2, alignment=0, spaceAfter=2)
    CAP = ParagraphStyle("CAP", parent=B, fontSize=8.5, leading=11, textColor=colors.HexColor("#555555"), alignment=0)
    TITLE = ParagraphStyle("T", parent=ss["Title"], textColor=colors.HexColor(NAVY), fontSize=22, leading=27, alignment=0)
    SUB = ParagraphStyle("S", parent=B, fontSize=11, textColor=colors.HexColor("#444444"), alignment=0)
    CELL = ParagraphStyle("C", parent=B, fontSize=8.6, leading=11, alignment=0, spaceAfter=0)

    P = lambda t, s=B: Paragraph(t, s)
    bullets = lambda items: [Paragraph(i, BL, bulletText="\u2022") for i in items]
    doc = SimpleDocTemplate(pdf_path, pagesize=A4, leftMargin=2 * cm, rightMargin=2 * cm,
                            topMargin=1.8 * cm, bottomMargin=1.8 * cm, title=REPORT_TITLE)
    W = A4[0] - 4 * cm

    def table(rows, widths):
        head = [Paragraph(f'<font color="white"><b>{esc(c)}</b></font>', CELL) for c in rows[0]]
        body = [[Paragraph(esc(c), CELL) for c in r] for r in rows[1:]]
        t = Table([head] + body, colWidths=widths, repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(NAVY)),
            ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#C9CED6")),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 3), ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F3F6FA")])]))
        return t

    def fig(path, w, cap):
        iw, ih = PILImage.open(path).size
        return KeepTogether([Image(path, width=w, height=w * ih / iw), P(cap, CAP), Spacer(1, 6)])

    def footer(c, d):
        c.saveState(); c.setFont("Helvetica", 8); c.setFillColor(colors.HexColor("#777777"))
        c.drawString(2 * cm, 1.1 * cm, REPORT_TITLE)
        c.drawRightString(A4[0] - 2 * cm, 1.1 * cm, f"Page {d.page}")
        c.restoreState()

    n_pat, n_mut = res["n_patients"], res["n_mutations"]
    n_non, n_cod = res["n_noncoding"], res["n_coding"]
    sp = res["spectrum"]; top_sub = sp["percent"].idxmax()
    rec = res["recurrent"]
    s = []

    # --- title page block
    s += [Spacer(1, 2 * cm), P(esc(REPORT_TITLE), TITLE), P(esc(REPORT_SUBTITLE), SUB), Spacer(1, 14)]
    kpi_head = ["MAF files", "Patients analysed", "Mutations (de-duplicated)", "Noncoding fraction"]
    kpi_val = [f"{info['n_files']:,}", f"{n_pat:,}", f"{n_mut:,}", f"{100 * n_non / n_mut:.1f}% ({n_non:,})"]
    t = Table([[Paragraph(f'<font size=8 color="#555555">{h}</font>', CELL) for h in kpi_head],
               [Paragraph(f'<font size=13 color="{NAVY}"><b>{esc(v)}</b></font>', CELL) for v in kpi_val]],
              colWidths=[W / 4] * 4)
    t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor(NAVY)),
                           ("INNERGRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#C9CED6")),
                           ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#F3F6FA")),
                           ("TOPPADDING", (0, 0), (-1, -1), 8), ("BOTTOMPADDING", (0, 0), (-1, -1), 8)]))
    s += [t, Spacer(1, 12)]

    # --- summary
    lead = rec.iloc[0] if len(rec) else None
    s += [P("Summary", H1),
          P("This report asks whether somatic mutations occur in melanoma, and which mutations are most frequent, "
            "as a first step toward a search for noncoding mutations. Mutation annotation format (MAF) files were "
            "merged, reduced to one tumour sample per patient, de-duplicated and classified by genomic region.")]
    items = [f"<b>Mutations are abundant:</b> {n_mut:,} somatic mutations in {n_pat:,} patients "
             f"(median {int(res['per_patient'].median()):,} per patient).",
             f"<b>Substitution spectrum:</b> {esc(top_sub)} accounts for {sp.loc[top_sub, 'percent']}% of coding point mutations."]
    if lead is not None:
        items.append(f"<b>Most recurrent protein change:</b> {esc(lead['Hugo_Symbol'])} {esc(lead['HGVSp_Short'])} in "
                     f"{lead['patients']} of {n_pat:,} patients ({lead['percent']}%).")
    items.append(f"<b>Noncoding mutations</b> make up {100 * n_non / n_mut:.1f}% of the data ({n_non:,}). Recurrent "
                 "noncoding positions are a preview only: no background mutation model has been applied.")
    s += bullets(items)

    # --- methods
    st = info["sample_types"]
    st_rows = [["Sample type", "Patients"]] + [[SAMPLE_NAME.get(k, f"Other ({k})"), f"{v:,}"] for k, v in st.items()]
    builds = ", ".join(sorted(res["builds"])) if res["builds"] else "not reported"
    s += [P("1. Data and methods", H1),
          P(f"<b>Input.</b> {info['n_files']:,} MAF files; {info['rows_read']:,} rows read. Genome build: {esc(builds)}."),
          P(f"<b>Sample selection.</b> The patient is the first 12 characters of the TCGA barcode. "
            f"{info['multi_sample_patients']} patient(s) had more than one sample in the files; one sample was kept "
            f"per patient, preferring primary (01), then metastatic (06). Duplicate mutation rows within a patient "
            f"were removed, leaving {n_mut:,} mutations in {n_pat:,} patients."),
          P("<b>Classification.</b> Variant_Classification was mapped to: coding (missense, silent, nonsense, "
            "frame-shift, in-frame, translation start site, nonstop), splice, intron, ncRNA (class 'RNA'), UTR, "
            "flank and intergenic."),
          Spacer(1, 4), table(st_rows, [W * 0.6, W * 0.4]), Spacer(1, 8)]
    vt = res["variant_types"]
    s += [P("Variant types", H2),
          table([["Type", "Count", "Share"]] + [[k, f"{v:,}", f"{100 * v / n_mut:.2f}%"] for k, v in vt.items()],
                [W * 0.4, W * 0.3, W * 0.3])]

    # --- section 2
    ct = res["class_table"]
    s += [PageBreak(), P("2. Where do mutations occur?", H1),
          P(f"{n_cod:,} of {n_mut:,} mutations ({100 * n_cod / n_mut:.1f}%) are coding; {n_non:,} "
            f"({100 * n_non / n_mut:.1f}%) are noncoding."),
          fig(F["region"], W * 0.9, "Figure 1. Mutations by region group (log scale)."),
          fig(F["share"], W * 0.7, "Figure 2. Coding versus noncoding share."),
          P("Variant classification", H2),
          table([["Classification", "Mutations", "Percent", "Region group"]] +
                [[r.classification, f"{r.mutations:,}", f"{r.percent:.2f}", r.region_group] for r in ct.itertuples()],
                [W * 0.34, W * 0.22, W * 0.18, W * 0.26]), Spacer(1, 6),
          fig(F["classes"], W * 0.95, "Figure 3. Top 12 variant classifications, coloured by region group."),
          fig(F["burden"], W * 0.85, "Figure 4. Distribution of mutations per patient (median marked).")]

    # --- section 3
    s += [PageBreak(), P("3. Which mutations are most frequent?", H1), P("3.1 Substitution spectrum", H2),
          table([["Substitution", "Count", "Percent"]] + [[i, f"{int(r['count']):,}", r["percent"]] for i, r in sp.iterrows()],
                [W * 0.4, W * 0.3, W * 0.3]), Spacer(1, 6),
          P(f"Substitutions are folded onto the pyrimidine strand (G&gt;A is counted as C&gt;T). In melanoma a large "
            f"C&gt;T share is the signature of UV damage."),
          fig(F["spectrum"], W * 0.85, "Figure 5. Substitution spectrum of coding point mutations."),
          P("3.2 Most recurrent protein changes", H2),
          table([["Gene", "Protein change", f"Patients (of {n_pat:,})", "Percent"]] +
                [[r.Hugo_Symbol, r.HGVSp_Short, r.patients, r.percent] for r in rec.itertuples()],
                [W * 0.22, W * 0.30, W * 0.28, W * 0.20]), Spacer(1, 6),
          fig(F["recurrent"], W * 0.95, "Figure 6. Most recurrent protein-level changes (BRAF and NRAS highlighted)."),
          KeepTogether([P("3.3 Genes mutated in the most patients", H2),
                        fig(F["genes"], W * 0.9, "Figure 7. Top genes by % of patients with a non-silent coding mutation."),
                        P("Very large genes (TTN, MUC16, ...) top such lists because they offer more positions to "
                          "mutate, and melanoma has a very high mutation rate. A frequently mutated gene is not "
                          "automatically a cancer driver; use a gene-length- and context-corrected test "
                          "(MutSigCV, dNdScv, OncodriveFML) to call drivers.")])]

    s += [P("3.4 Driver hotspots examined in detail", H2)]
    for g in DRIVER_GENES:
        tab = res["gene_tables"].get(g)
        n_g = res["gene_patients"].get(g, 0)
        if tab is None or tab.empty:
            s += [P(f"<b>{esc(g)}</b>: no coding mutations with protein positions found.")]
            continue
        s += [P(f"<b>{esc(g)}</b>: {n_g} patients ({100 * n_g / n_pat:.1f}%) have a non-silent coding mutation."),
              table([["Protein position", "Classification", "Patients"]] +
                    [[r.pos, r.Variant_Classification, r.patients] for r in tab.itertuples()],
                    [W * 0.3, W * 0.4, W * 0.3]), Spacer(1, 6)]
        if g == "BRAF" and str(tab.iloc[0]["pos"]) == "640":
            s += [P("Note: position 640 is how this annotation numbers the BRAF hotspot on its longer transcript; it "
                    "corresponds to the widely quoted V600 hotspot.", CAP)]

    # --- section 4
    s += [P("4. Noncoding mutations (preview)", H1)]
    rc = res["region_counts"]
    nc_parts = ", ".join(f"{k} {v:,}" for k, v in rc.items() if k != "coding")
    s += [P(f"Noncoding categories hold {n_non:,} mutations ({100 * n_non / n_mut:.1f}%): {esc(nc_parts)}. "
            "The following is exploratory only.")]
    if "ncrna" in F:
        s += [KeepTogether([P("4.1 ncRNA-class genes", H2),
                            fig(F["ncrna"], W * 0.85, "Figure 8. Genes with the most mutations in the 'RNA' class."),
                            P("Pseudogenes and repetitive loci are prone to mapping artefacts, so these genes need "
                              "filtering before any functional interpretation.")])]
    nr = res["nc_recurrent"]
    s += [P(f"4.2 Noncoding positions mutated in {MIN_PATIENTS_NONCODING} or more patients", H2)]
    if nr.empty:
        s += [P("No noncoding position reached the recurrence cut-off.")]
    else:
        s += [table([["Position (genome build as in files)", "Class", "Gene", "Patients"]] +
                    [[r.position, r.Variant_Classification, r.Hugo_Symbol, r.patients] for r in nr.itertuples()],
                    [W * 0.34, W * 0.24, W * 0.28, W * 0.14])]
    s += [Spacer(1, 6),
          P("<b>Caution.</b> Melanoma has a very high UV mutation rate, so many positions recur by chance. Real "
            "hotspot calling needs a background mutation model (trinucleotide context, replication timing, local "
            "mutation rate) and multiple-testing correction.")]

    # --- limitations / next steps
    s += [P("5. Limitations", H1)] + bullets([
        "<b>Sequencing scope.</b> If most mutations are coding and noncoding calls sit only at splice, UTR and flank "
        "positions near genes, the data are probably whole-exome. Deep intronic and intergenic regions would then be "
        "unsampled. Confirm the experimental strategy (WXS or WGS) in the GDC manifest.",
        "<b>Mixed sample types.</b> Primary and metastatic samples are pooled (see the table in section 1).",
        "<b>Open-access (masked) files.</b> Some variants may be masked or filtered compared with controlled-access data.",
        "<b>No statistical test.</b> Frequency lists reflect gene size and mutation rate; no driver significance is computed.",
        "<b>Annotation transcript.</b> Protein positions depend on the transcript used by the annotation."])
    s += [P("6. Suggested next steps", H1)] + bullets([
        "Confirm WXS versus WGS; use WGS if true noncoding discovery is the goal.",
        "Run a driver test for coding genes (MutSigCV or dNdScv).",
        "Build a background model for noncoding hotspots and correct for multiple testing.",
        "Filter pseudogenes and low-mappability regions; compare with known noncoding events (for example TERT promoter).",
        "Stratify by primary versus metastatic samples."])

    doc.build(s, onFirstPage=footer, onLaterPages=footer)


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maf-dir", default=DEFAULT_MAF_DIR, help="folder containing the MAF files")
    ap.add_argument("--out-dir", default=DEFAULT_OUT_DIR, help="folder for the PDF, figures and tables")
    args = ap.parse_args()

    files = find_maf_files(args.maf_dir)
    if not files:
        sys.exit(f"No MAF files found under: {args.maf_dir}")
    print(f"MAF files found: {len(files)}")
    fig_dir, tab_dir = os.path.join(args.out_dir, "figures"), os.path.join(args.out_dir, "tables")
    os.makedirs(fig_dir, exist_ok=True); os.makedirs(tab_dir, exist_ok=True)

    raw = load_mafs(files)
    df, info = clean(raw)
    info["n_files"] = len(files)
    del raw
    print(f"Patients: {info['n_patients']}   Mutations: {info['n_mutations']:,}")

    coding = df[df["region_group"] == "coding"]
    nonsilent = coding[coding["Variant_Classification"] != "Silent"]
    noncoding = df[df["region_group"] != "coding"]
    n_pat = info["n_patients"]

    class_table = (df.groupby(["Variant_Classification", "region_group"]).size()
                     .reset_index(name="mutations").sort_values("mutations", ascending=False))
    class_table["percent"] = 100 * class_table["mutations"] / len(df)
    class_table = class_table.rename(columns={"Variant_Classification": "classification"})

    ncrna = df[df["Variant_Classification"] == "RNA"]
    ncrna_genes = (ncrna.groupby("Hugo_Symbol").size().sort_values(ascending=False).head(15)
                        .reset_index(name="mutations"))

    res = {
        "n_patients": n_pat, "n_mutations": len(df), "n_coding": len(coding), "n_noncoding": len(noncoding),
        "region_counts": df["region_group"].value_counts(),
        "class_table": class_table,
        "variant_types": df["Variant_Type"].value_counts() if "Variant_Type" in df else pd.Series(dtype=int),
        "builds": set(df["NCBI_Build"].dropna().unique()) if "NCBI_Build" in df else set(),
        "spectrum": pyrimidine_spectrum(coding),
        "recurrent": top_recurrent_protein(nonsilent, n_pat),
        "top_genes": top_genes(nonsilent, n_pat),
        "ncrna_genes": ncrna_genes,
        "nc_recurrent": noncoding_recurrent(noncoding, MIN_PATIENTS_NONCODING),
        "per_patient": df.groupby("patient").size(),
        "gene_tables": {g: gene_positions(coding, g) for g in DRIVER_GENES},
        "gene_patients": {g: nonsilent[nonsilent["Hugo_Symbol"] == g]["patient"].nunique() for g in DRIVER_GENES},
    }

    # save every table as CSV
    res["class_table"].to_csv(os.path.join(tab_dir, "variant_classes.csv"), index=False)
    res["region_counts"].to_csv(os.path.join(tab_dir, "region_groups.csv"), header=["mutations"])
    res["spectrum"].to_csv(os.path.join(tab_dir, "substitution_spectrum.csv"))
    res["recurrent"].to_csv(os.path.join(tab_dir, "recurrent_protein_changes.csv"), index=False)
    res["top_genes"].to_csv(os.path.join(tab_dir, "top_genes.csv"), index=False)
    res["nc_recurrent"].to_csv(os.path.join(tab_dir, "noncoding_recurrent_positions.csv"), index=False)
    res["per_patient"].to_csv(os.path.join(tab_dir, "mutations_per_patient.csv"), header=["mutations"])
    for g, t in res["gene_tables"].items():
        t.to_csv(os.path.join(tab_dir, f"{g}_positions.csv"), index=False)

    figs = make_figures(df, res, fig_dir)
    pdf = os.path.join(args.out_dir, "Mutation_Landscape_Report.pdf")
    build_pdf(pdf, res, info, figs)
    print(f"\nDone. Report: {pdf}")


if __name__ == "__main__":
    main()