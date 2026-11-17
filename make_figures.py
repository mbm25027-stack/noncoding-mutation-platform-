#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_figures.py -- Figures 1-5 for the NC-MutScan manuscript
=============================================================

Run from the project folder (Windows PowerShell, venv active):

    python make_figures.py                    # all figures, PNG + PDF, 300 dpi
    python make_figures.py --only 2 3         # selected figures
    python make_figures.py --formats png pdf svg tiff --outdir paper_figures
    python make_figures.py --dpi 600          # higher raster resolution

Requires only: numpy, matplotlib.

How data is handled
-------------------
* Every number that appears in the experiment reports (EXP 01-10) is embedded
  below, so every panel can be drawn without any input file.
* Four panels are *better* with per-variant data. If the optional CSV files
  listed in OPTIONAL_DATA exist, the script uses them (real ROC curves, the
  per-variant false-negative scatter, the TERT per-position track). If they do
  not, the script draws a clearly labelled summary version of that panel from
  the report numbers and prints a notice, so nothing is silently faked.
* Built-in self-checks recompute AUC (when scores are supplied), post-test PPV
  and Benjamini-Hochberg q-values and compare them with the values in the
  reports. A mismatch prints a WARNING.

Output goes to ./paper_figures (never overwrites your existing Figure*.png).

Style notes (v2)
----------------
Only presentation was changed relative to v1; no number, statistic, threshold
or label of scientific meaning was altered. Raster output is 300 dpi by
default; vector formats (PDF/SVG) keep editable, embedded text.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
from matplotlib.gridspec import GridSpec                           # noqa: E402
from matplotlib.patches import Rectangle                           # noqa: E402

HERE = Path(__file__).resolve().parent

# --------------------------------------------------------------------------
# OPTIONAL per-variant inputs (all optional; see the schema in each comment)
# --------------------------------------------------------------------------
OPTIONAL_DATA = {
    # Fig 3a. n = 139 blind variants with SpliceAI. One row per variant.
    # Columns (higher score = more pathogenic for every score column):
    #   label (1 pathogenic / 0 benign), spliceai, cadd_phred, maxentscan,
    #   exon_junction (use -|d_junction|), nc_s_raw, nc_s_unified
    "blind_scores": HERE / "figdata" / "blind_scores.csv",
    # Fig 4b. The 38 false negatives of the blind set.
    #   Columns: cls (1-7), wt_score, dpwm_splice, spliceai_ge_020 (0/1)
    "fn_variants": HERE / "figdata" / "fn_variants.csv",
    # Fig 5b. All 753 TERT core-promoter SNVs (GRCh38 chr5:1,294,995-1,295,245).
    #   Columns: pos, mpra_effect (log2), activating (0/1), call_denovo (0/1)
    "tert_snvs": HERE / "figdata" / "tert_snvs.csv",
}
SCORE_COLUMNS = {      # key -> (CSV column, legend label, reported AUC)
    "spliceai":      ("spliceai",       "SpliceAI",                0.946),
    "cadd":          ("cadd_phred",     "CADD v1.7",               0.871),
    "maxent":        ("maxentscan",     "MaxEntScan",              0.832),
    "junction":      ("exon_junction",  "Exon-junction proximity", 0.819),
    "nc_raw":        ("nc_s_raw",       "NC-MutScan (S_raw)",      0.785),
    "nc_unified":    ("nc_s_unified",   "NC-MutScan (S_unified)",  0.733),
}

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
MM = 1 / 25.4
W1, W2 = 89 * MM, 183 * MM                       # Nature single / double column
# Okabe-Ito colour-blind-safe palette (unchanged hues, used sparingly)
C = dict(blue="#0072B2", orange="#E69F00", green="#009E73", red="#D55E00",
         sky="#56B4E9", pink="#CC79A7", grey="#7F7F7F", dark="#1A1A1A",
         light="#EDEDED", mid="#BDBDBD", fill="#F5F5F5", text2="#404040")
FS_TICK, FS_NOTE = 6.0, 6.0          # nothing smaller than 6 pt anywhere


def set_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
        "axes.titleweight": "normal", "axes.labelpad": 2.5,
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK, "legend.fontsize": FS_TICK,
        "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.pad": 2, "ytick.major.pad": 2,
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.0, "lines.markeredgewidth": 0.5,
        "patch.linewidth": 0.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.dpi": 300, "figure.dpi": 150, "mathtext.default": "regular",
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "legend.frameon": False, "legend.handletextpad": 0.5,
    })


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """Wilson score interval for a proportion (returned as fractions)."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def ppv_from_counts(tp: int, npos: int, tn: int, nneg: int, prior):
    sens, spec = tp / npos, tn / nneg
    prior = np.asarray(prior, dtype=float)
    return sens * prior / (sens * prior + (1 - spec) * (1 - prior))


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / (np.arange(m) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(m)
    q[order] = np.minimum(ranked, 1.0)
    return q


def panel_label(ax, text: str, dx: float = -0.16, dy: float = 1.06) -> None:
    ax.text(dx, dy, text, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left")


def fig_label(fig, text: str, x: float, y: float) -> None:
    """Bold lowercase panel letter at figure coordinates (0-1)."""
    fig.text(x, y, text, fontsize=8, fontweight="bold", va="top", ha="left")


def save(fig, name: str, outdir: Path, formats) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = outdir / f"{name}.{fmt}"
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if fmt == "tiff" else {}
        fig.savefig(path, format=fmt, bbox_inches=None, **kw)
        print(f"  wrote {path}")
    plt.close(fig)


def read_csv(path: Path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def warn(msg: str) -> None:
    print(f"  WARNING: {msg}")


def note(msg: str) -> None:
    print(f"  note: {msg}")


def roc_curve(y: np.ndarray, s: np.ndarray):
    """ROC points (fpr, tpr) and AUC; ties handled by distinct thresholds."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, float)
    order = np.argsort(-s, kind="mergesort")
    y, s = y[order], s[order]
    distinct = np.where(np.diff(s))[0]
    idx = np.r_[distinct, len(y) - 1]
    tps = np.cumsum(y)[idx]
    fps = 1 + idx - tps
    tpr = np.r_[0, tps / y.sum()]
    fpr = np.r_[0, fps / (len(y) - y.sum())]
    auc = float(np.trapezoid(tpr, fpr)) if hasattr(np, "trapezoid") else float(np.trapz(tpr, fpr))
    return fpr, tpr, auc


def bar_with_ci(ax, xs, vals, lo, hi, width, color, label=None):
    """Bars with Wilson-CI whiskers drawn on top of the bar."""
    ax.bar(xs, vals, width, color=color, label=label, lw=0, zorder=2)
    ax.errorbar(xs, vals, yerr=[lo, hi], fmt="none", ecolor=C["dark"],
                elinewidth=0.6, capsize=1.5, capthick=0.6, zorder=3)


# --------------------------------------------------------------------------
# Numbers from the experiment reports
# --------------------------------------------------------------------------
# EXP 04 -- runtime (microseconds per variant) and 54-locus concordance
RUNTIME = [("Window preparation", 58.2), ("Splice + polyA modules", 164.5),
           ("TF scan + AID/SHM masking", 273.0)]
LOCI = [("Canonical hotspots", 4), ("Endogenous AID/SHM regions", 10),
        ("Intergenic controls", 36), ("Truncated inputs", 4)]

# EXP 05 -- locked blind set (chr12-22), counts
BLIND = dict(tp=37, fp=1, tn=74, fn=38)
ABLATION = [("Full pipeline", 37, 1, 74, 38),
            ("Splice/polyA only", 36, 0, 75, 39),
            ("TF motifs only", 1, 1, 74, 74)]
# EXP 06 -- LOCO (N = 300; 150 P / 150 B)
LOCO_FIXED = dict(tp=67, fn=83, tn=141, fp=9)
LOCO_NESTED = dict(tp=55, fn=95, tn=138, fp=12)
# EXP 06 -- matched-specificity sensitivity (69 pathogenic, ~98.6% specificity)
MATCHED = [("SpliceAI", 53), ("MaxEntScan", 36), ("Exon-junction proximity", 34),
           ("NC-MutScan", 33), ("CADD v1.7", 27)]
AUC_REPORTED = [(v[1], v[2]) for v in SCORE_COLUMNS.values()]
# EXP 07 -- false-negative classes: (class, label, n, SpliceAI-detected, median WT, median |dPWM|)
FN_CLASSES = [
    (1, "Extended splice, sub-threshold shift", 15, 13, 6.45, 1.42),
    (2, "Extended splice, outside PWM window", 7, 5, 0.00, 0.00),
    (3, "Canonical \u00b11\u20132 bp, weak wild-type", 11, 11, 4.18, 4.61),
    (4, "Deep-intronic / distal regulatory", 2, 0, None, 1.10),      # WT < 5.0
    (5, "Deep-intronic pseudoexon", 1, 1, 3.10, 2.95),
    (6, "Putative branch point", 1, 0, 0.00, 0.00),
    (7, "AID/SHM-masked splice variant", 1, 1, 6.80, 3.40),
]
# EXP 08 -- TERT core promoter, 4 activating / 749 neutral SNVs
TERT_MODES = [("High-stringency\n(ClinVar cut)", 0, 24, 725, 4),
              ("De novo TF-gain\nmode", 1, 26, 723, 3)]
# EXP 09 -- 29 Kircher libraries.
# (name, role, TP, calls, LR+, lo, hi, Fisher P, Haldane-Anscombe corrected)
LIBS = [
    ("BCL11A", "ind", 1, 87, 1.86, 0.29, 11.74, 0.430, False),
    ("F9", "ind", 5, 18, 2.41, 0.89, 6.53, 0.086, False),
    ("FOXE1", "ind", 0, 11, 27.26, 2.31, 321.79, 1.000, True),
    ("GP1BA", "ind", 3, 27, 0.70, 0.21, 2.28, 0.785, False),
    ("HBB", "ind", 0, 15, 0.50, 0.03, 7.95, 1.000, True),
    ("HBG1", "ind", 2, 3, 5.56, 0.51, 60.75, 0.172, False),
    ("HNF4A", "ind", 0, 18, 1.52, 0.10, 23.09, 1.000, True),
    ("IRF4", "ind", 10, 24, 0.57, 0.26, 1.27, 0.210, False),
    ("IRF6", "ind", 21, 74, 1.35, 0.83, 2.18, 0.246, False),
    ("LDLR", "ind", 6, 11, 1.49, 0.46, 4.84, 0.551, False),
    ("LDLR.2", "ind", 6, 12, 1.14, 0.37, 3.51, 1.000, False),
    ("MSMB", "ind", 8, 62, 4.04, 2.11, 7.74, 8.62e-4, False),
    ("MYC rs11986220", "ind", 0, 67, 2.07, 0.16, 26.28, 1.000, True),
    ("MYC rs6983267", "ind", 3, 57, 2.41, 0.80, 7.21, 0.134, False),
    ("PKLR-24h", "ind", 5, 23, 0.58, 0.22, 1.55, 0.365, False),
    ("PKLR-48h", "ind", 5, 24, 0.46, 0.17, 1.23, 0.133, False),
    ("RET", "ind", 2, 40, 2.53, 0.65, 9.81, 0.194, False),
    ("SORT1", "ind", 15, 39, 0.52, 0.27, 0.97, 0.048, False),
    ("SORT1-flip", "ind", 10, 42, 0.41, 0.20, 0.82, 0.010, False),
    ("SORT1.2", "ind", 14, 44, 0.52, 0.28, 0.96, 0.044, False),
    ("TCF7L2", "ind", 1, 41, 0.83, 0.12, 5.88, 1.000, False),
    ("TERT-GAa", "tert", 8, 15, 2.19, 0.81, 5.89, 0.162, False),
    ("TERT-GBM", "tert", 12, 18, 2.10, 0.80, 5.49, 0.150, False),
    ("TERT-GSc", "tert", 11, 17, 2.50, 0.94, 6.62, 0.077, False),
    ("TERT-HEK", "tert", 3, 18, 0.88, 0.26, 2.95, 1.000, False),
    ("UC88", "ind", 3, 55, 2.32, 0.77, 6.98, 0.144, False),
    ("ZFAND3", "ind", 3, 59, 0.33, 0.10, 1.04, 0.051, False),
    ("ZRSh-13", "ind", 0, 55, 1.59, 0.11, 22.26, 1.000, True),
    ("ZRSh-13h2", "ind", 0, 51, 0.83, 0.06, 12.61, 1.000, True),
]
POOLED = [("All 29 libraries", 157, 1027, 0.78, 0.66, 0.92, 2.50e-3),
          ("25 non-TERT libraries", 123, 959, 0.69, 0.58, 0.83, 5.18e-5)]


# ==========================================================================
# FIGURE 1  --  pipeline, runtime, 54-locus verification
# ==========================================================================
def fig1(outdir, formats):
    fig = plt.figure(figsize=(W2, 112 * MM))
    axa = fig.add_axes([0.02, 0.46, 0.96, 0.51])
    axb = fig.add_axes([0.05, 0.14, 0.45, 0.19])
    axc = fig.add_axes([0.745, 0.14, 0.235, 0.21])

    # ---- a: schematic ----------------------------------------------------
    axa.set_xlim(0, 100)
    axa.set_ylim(0, 44)
    axa.axis("off")
    fig_label(fig, "a", 0.004, 0.985)

    def box(x, y, w, h, title, sub="", fc="#FFFFFF", accent=None, tfs=6.8, sfs=6.0):
        """Flat rectangular box; optional thin accent bar on the left edge."""
        axa.add_patch(Rectangle((x, y), w, h, fc=fc, ec=C["dark"], lw=0.6, zorder=2))
        if accent:
            axa.add_patch(Rectangle((x, y), 1.1, h, fc=accent, ec="none", zorder=3))
        tx = x + w / 2 + (0.5 if accent else 0.0)
        if sub:
            axa.text(tx, y + h * 0.72, title, ha="center", va="center", fontsize=tfs,
                     fontweight="bold", color=C["dark"], zorder=4)
            axa.text(tx, y + h * 0.31, sub, ha="center", va="center", fontsize=sfs,
                     color=C["text2"], linespacing=1.3, zorder=4)
        else:
            axa.text(tx, y + h / 2, title, ha="center", va="center", fontsize=tfs,
                     fontweight="bold", color=C["dark"], linespacing=1.3, zorder=4)

    def arrow(x0, y0, x1, y1):
        axa.annotate("", xy=(x1, y1), xytext=(x0, y0),
                     arrowprops=dict(arrowstyle="-|>", lw=0.6, color=C["dark"],
                                     shrinkA=0, shrinkB=0, mutation_scale=6))

    mid = 20.5
    axa.text(50, 43.0, "Variant-centred WT vs mutant differential scoring (\u0394PWM = MUT \u2212 WT, in bits)",
             ha="center", va="center", fontsize=6.8, color=C["dark"])
    box(1, mid - 5.5, 17, 11, "Wild-type and\nmutant sequence", fc=C["fill"])
    axa.text(9.5, mid - 7.6, "\u00b125 bp window overlapping\nthe mutated base, both strands",
             ha="center", va="top", fontsize=6.0, color=C["text2"], linespacing=1.25)
    mods = [
        (28.0, "Splice junctions",
         "U1 donor 9-mer \u00b7 U2 acceptor + polypyrimidine tract\n|\u0394PWM| \u2265 2.8 bits; WT \u2265 5.6 (5\u2032) / 6.4 (3\u2032)", C["blue"]),
        (15.0, "3\u2032-UTR polyadenylation",
         "AATAAA / ATTAAA loss in \u2265 65% A+T window\n4.2-bit floor", C["orange"]),
        (2.0, "Transcription factors (9 PWMs)",
         "Sp1 GABPA MYB ETS1 E2F1 AP-1 GATA3 FOXA1 CTCF\n|\u0394PWM| \u2265 4.2 bits; active allele \u2265 7.8 bits", C["green"]),
    ]
    for y, t, s, acc in mods:
        box(26, y, 38, 11, t, s, fc="#FFFFFF", accent=acc)
        arrow(18.2, mid, 25.6, y + 5.5)
        arrow(64.4, y + 5.5, 69.4, mid)
    box(70, mid - 7, 15, 14, "AID/SHM mask", "WRCY / RGYW,\n\u2265 6 dual-strand hits,\nno dominant motif",
        fc="#FFFFFF", accent=C["red"], sfs=5.9)
    box(89, mid - 5.5, 10.5, 11, "S_unified", "\u2265 4.2 bits =\nrule-in call", fc=C["fill"], sfs=5.9)
    arrow(85.4, mid, 88.6, mid)
    axa.text(77.5, mid - 9.4, "Promoter de novo mode:\n\u0394PWM \u2265 +3.5, mutant \u2265 7.5 bits",
             ha="center", va="top", fontsize=6.0, color=C["text2"], linespacing=1.25)

    # ---- b: runtime -------------------------------------------------------
    total = sum(v for _, v in RUNTIME)
    cols = [C["sky"], C["blue"], C["green"]]
    short = ["Window\npreparation", "Splice + polyA\nmodules", "TF scan + AID/SHM\nmasking"]
    left = 0.0
    for (lab, v), col, sl in zip(RUNTIME, cols, short):
        axb.barh(0, v, left=left, height=0.55, color=col, ec="white", lw=0.8)
        axb.text(left + v / 2, 0, f"{v:.1f}", ha="center", va="center", fontsize=6.5,
                 color="white", fontweight="bold")
        axb.text(left + v / 2, -0.42, sl, ha="center", va="top", fontsize=6.0,
                 color=C["text2"], linespacing=1.15)
        left += v
    axb.set_xlim(0, total * 1.02)
    axb.set_ylim(-1.45, 1.0)
    axb.set_yticks([])
    axb.spines["left"].set_visible(False)
    axb.set_xlabel("CPU time per variant (\u00b5s, single thread)")
    axb.text(0, 0.5, f"Total {total:.1f} \u00b5s per variant \u2248 2,017 variants s$^{{-1}}$;  HTTP API < 15 ms",
             fontsize=6.3, ha="left", va="bottom", color=C["dark"])
    fig_label(fig, "b", 0.004, 0.40)

    # ---- c: 54 loci -------------------------------------------------------
    names = ["Canonical hotspots", "AID/SHM regions", "Intergenic controls", "Truncated inputs"][::-1]
    ns = [k for _, k in LOCI][::-1]
    axc.barh(names, ns, color=C["blue"], height=0.6, lw=0)
    for i, k in enumerate(ns):
        axc.text(k + 0.8, i, f"{k}/{k}", va="center", fontsize=6.2)
    axc.set_xlim(0, 46)
    axc.set_xlabel("GRCh38 loci (n; all concordant)")
    axc.set_title("54/54 concordant, 0 false calls in controls", fontsize=6.5, loc="right", pad=4)
    fig_label(fig, "c", 0.56, 0.40)
    save(fig, "Fig1_pipeline_runtime_verification", outdir, formats)


# ==========================================================================
# FIGURE 2  --  locked blind ClinVar benchmark
# ==========================================================================
def fig2(outdir, formats):
    fig = plt.figure(figsize=(W2, 66 * MM))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[0.8, 1.0, 1.15],
                  left=0.10, right=0.99, top=0.88, bottom=0.2, wspace=0.5)
    axa, axb, axc = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                     fig.add_subplot(gs[0, 2]))

    # ---- a: confusion matrix ----------------------------------------------
    b = BLIND
    # correct calls in colour, errors in neutral grey (same layout as before)
    cells = [(0, 1, b["tp"], C["blue"], "TP"), (1, 1, b["fn"], C["light"], "FN"),
             (0, 0, b["fp"], C["light"], "FP"), (1, 0, b["tn"], C["blue"], "TN")]
    # grid: x = predicted (0 rule-in, 1 not called); y = truth (1 pathogenic top, 0 benign bottom)
    for xi, yi, v, col, tag in cells:
        axa.add_patch(Rectangle((xi, yi), 1, 1, fc=col, ec="white", lw=1.5))
        strong = col == C["blue"]
        axa.text(xi + 0.5, yi + 0.58, str(v), ha="center", va="center", fontsize=12,
                 fontweight="bold", color="white" if strong else C["dark"])
        axa.text(xi + 0.5, yi + 0.25, tag, ha="center", va="center", fontsize=6.5,
                 color="white" if strong else C["dark"])
    axa.set_xlim(0, 2)
    axa.set_ylim(0, 2)
    axa.set_xticks([0.5, 1.5])
    axa.set_xticklabels(["Rule-in", "Not called"])
    axa.set_yticks([0.5, 1.5])
    axa.set_yticklabels(["Benign\n(n = 75)", "Pathogenic\n(n = 75)"])
    axa.set_xlabel("NC-MutScan call")
    axa.tick_params(length=0)
    for s in axa.spines.values():
        s.set_visible(False)
    axa.set_title("Locked blind set, chr12\u201322", fontsize=6.8, pad=4)
    fig_label(fig, "a", 0.004, 0.985)

    # ---- b: post-test PPV vs prior ----------------------------------------
    prior = np.logspace(-3, -1, 200)
    curves = [("Blind set (n = 150)", (37, 75, 74, 75), C["red"]),
              ("LOCO, fixed (n = 300)", (67, 150, 141, 150), C["blue"]),
              ("LOCO, nested tuning", (55, 150, 138, 150), C["grey"])]
    for lab, (tp, npos, tn, nneg), col in curves:
        axb.plot(prior * 100, ppv_from_counts(tp, npos, tn, nneg, prior) * 100,
                 color=col, lw=1.2, label=lab)
    pts = np.array([0.001, 0.005, 0.01, 0.02, 0.05, 0.10])
    axb.plot(pts * 100, ppv_from_counts(37, 75, 74, 75, pts) * 100, "o", ms=3.0,
             mfc=C["red"], mec="white", mew=0.4, ls="none")
    p3 = np.array([0.01, 0.05, 0.10])
    axb.plot(p3 * 100, ppv_from_counts(67, 150, 141, 150, p3) * 100, "o", ms=3.0,
             mfc=C["blue"], mec="white", mew=0.4, ls="none")
    # self-check against report values
    rep_blind = [3.57, 15.66, 27.21, 43.02, 66.07, 80.43]
    got = ppv_from_counts(37, 75, 74, 75, pts) * 100
    if np.max(np.abs(got - rep_blind)) > 0.05:
        warn(f"blind PPV differs from report: {np.round(got, 2)} vs {rep_blind}")
    rep_loco = [7.0, 28.1, 45.3]
    got = ppv_from_counts(67, 150, 141, 150, p3) * 100
    if np.max(np.abs(got - rep_loco)) > 0.15:
        warn(f"LOCO PPV differs from report: {np.round(got, 2)} vs {rep_loco}")
    axb.set_xscale("log")
    axb.set_xlim(0.1, 10)
    axb.set_ylim(0, 100)
    axb.set_xticks([0.1, 0.5, 1, 2, 5, 10])
    axb.set_xticklabels(["0.1", "0.5", "1", "2", "5", "10"])
    axb.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    axb.set_xlabel("Prior probability of pathogenicity (%)")
    axb.set_ylabel("Post-test PPV (%)")
    axb.legend(loc="upper left", handlelength=1.4)
    axb.set_title("Prior-adjusted PPV", fontsize=6.8, pad=4)
    fig_label(fig, "b", 0.345, 0.985)

    # ---- c: ablation --------------------------------------------------------
    x = np.arange(len(ABLATION))
    w = 0.36
    for j, (metric, col) in enumerate([("Sensitivity", C["blue"]), ("Specificity", C["orange"])]):
        vals, lo, hi = [], [], []
        for _, tp, fp, tn, fn in ABLATION:
            k, n = (tp, tp + fn) if j == 0 else (tn, tn + fp)
            lo_, hi_ = wilson(k, n)
            vals.append(100 * k / n)
            lo.append(100 * k / n - 100 * lo_)
            hi.append(100 * hi_ - 100 * k / n)
        xs = x + (j - 0.5) * w
        bar_with_ci(axc, xs, vals, lo, hi, w, col, metric)
        for xx, v, h in zip(xs, vals, hi):
            axc.text(xx, v + h + 1.5, f"{v:.1f}", ha="center", va="bottom", fontsize=6.0)
    axc.set_xticks(x)
    axc.set_xticklabels([a[0].replace(" only", "\nonly").replace("Full ", "Full\n") for a in ABLATION])
    axc.set_ylim(0, 128)
    axc.set_yticks([0, 25, 50, 75, 100])
    axc.set_ylabel("Percent (95% Wilson CI)")
    axc.legend(loc="upper right", ncol=2, handlelength=1.0, columnspacing=0.8)
    axc.set_title("Component ablation", fontsize=6.8, pad=4)
    fig_label(fig, "c", 0.665, 0.985)
    save(fig, "Fig2_blind_clinvar_benchmark", outdir, formats)


# ==========================================================================
# FIGURE 3  --  comparison with established predictors, LOCO
# ==========================================================================
def fig3(outdir, formats):
    fig = plt.figure(figsize=(W2, 68 * MM))
    axa = fig.add_axes([0.20, 0.19, 0.19, 0.69])
    axb = fig.add_axes([0.555, 0.19, 0.14, 0.69])
    axc = fig.add_axes([0.86, 0.19, 0.13, 0.69])
    fig_label(fig, "a", 0.004, 0.985)
    fig_label(fig, "b", 0.415, 0.985)
    fig_label(fig, "c", 0.775, 0.985)
    colors = {"spliceai": C["blue"], "cadd": C["sky"], "maxent": C["green"],
              "junction": C["grey"], "nc_raw": C["red"], "nc_unified": C["orange"]}

    # ---- a: ROC (per-variant scores) or AUC summary ------------------------
    path = OPTIONAL_DATA["blind_scores"]
    drawn = False
    if path.exists():
        rows = read_csv(path)
        y = np.array([int(float(r["label"])) for r in rows])
        axa.plot([0, 1], [0, 1], ls="--", lw=0.6, color=C["mid"])
        for key, (col, lab, auc_rep) in SCORE_COLUMNS.items():
            if col not in rows[0]:
                warn(f"column '{col}' missing in {path.name}; skipping {lab}")
                continue
            s = np.array([float(r[col]) for r in rows])
            fpr, tpr, auc = roc_curve(y, s)
            if abs(auc - auc_rep) > 0.01:
                warn(f"{lab}: AUC {auc:.3f} from CSV vs {auc_rep:.3f} in report")
            axa.plot(fpr, tpr, lw=1.2 if "nc_" in key else 0.9, color=colors[key],
                     label=f"{lab} ({auc:.3f})")
        axa.axvline(1 - 0.9857, ls=":", lw=0.6, color=C["dark"])
        axa.set_xlim(0, 1)
        axa.set_ylim(0, 1.02)
        axa.set_aspect("equal", adjustable="box")
        axa.set_xlabel("1 \u2212 specificity")
        axa.set_ylabel("Sensitivity")
        axa.legend(loc="lower right", fontsize=6.0, handlelength=1.3, labelspacing=0.3)
        drawn = True
        n_note = f"n = {len(y)}"
    if not drawn:
        note(f"{path.name} not found; Fig 3a drawn as AUC summary (reported values). "
             "Add figdata/blind_scores.csv for real ROC curves.")
        keys = list(SCORE_COLUMNS)
        aucs = [SCORE_COLUMNS[k][2] for k in keys]
        order = np.argsort(aucs)
        for i, j in enumerate(order):
            k = keys[j]
            axa.plot([0.5, aucs[j]], [i, i], color=C["mid"], lw=0.8)
            axa.plot(aucs[j], i, "o", color=colors[k], ms=4.5, mec="white", mew=0.4)
            axa.text(aucs[j] + 0.016, i, f"{aucs[j]:.3f}", va="center", fontsize=6.0)
        axa.set_yticks(range(len(keys)))
        axa.set_yticklabels([SCORE_COLUMNS[keys[j]][1] for j in order])
        axa.set_xlim(0.5, 1.02)
        axa.set_xlabel("ROC-AUC")
        n_note = "n = 139"
    axa.set_title(f"Discrimination ({n_note})", fontsize=6.8, pad=4)

    # ---- b: sensitivity at matched specificity -----------------------------
    order = sorted(MATCHED, key=lambda t: t[1])
    for i, (lab, k) in enumerate(order):
        lo, hi = wilson(k, 69)
        col = C["red"] if lab == "NC-MutScan" else C["blue"]
        axb.barh(i, 100 * k / 69, color=col, height=0.62, lw=0, zorder=2)
        axb.errorbar(100 * k / 69, i, xerr=[[100 * k / 69 - 100 * lo], [100 * hi - 100 * k / 69]],
                     fmt="none", ecolor=C["dark"], elinewidth=0.6, capsize=1.5, capthick=0.6, zorder=3)
        axb.text(101, i, f"{100 * k / 69:.1f}% ({k}/69)", va="center", fontsize=6.0)
    axb.set_yticks(range(len(order)))
    axb.set_yticklabels([o[0] for o in order])
    axb.set_xlim(0, 100)
    axb.set_xlabel("Sensitivity at ~98.6% specificity (%)")
    axb.set_title("Matched operating point", fontsize=6.8, pad=4)

    # ---- c: LOCO fixed vs nested --------------------------------------------
    labs = ["Fixed\nbuild", "Nested\n\u03c4_TF tuning"]
    x = np.arange(2)
    w = 0.36
    for j, (metric, col) in enumerate([("Sensitivity", C["blue"]), ("Specificity", C["orange"])]):
        vals, lo, hi = [], [], []
        for d in (LOCO_FIXED, LOCO_NESTED):
            k, n = (d["tp"], d["tp"] + d["fn"]) if j == 0 else (d["tn"], d["tn"] + d["fp"])
            lo_, hi_ = wilson(k, n)
            vals.append(100 * k / n)
            lo.append(100 * k / n - 100 * lo_)
            hi.append(100 * hi_ - 100 * k / n)
        xs = x + (j - 0.5) * w
        bar_with_ci(axc, xs, vals, lo, hi, w, col, metric)
        for xx, v, h in zip(xs, vals, hi):
            axc.text(xx, v + h + 1.5, f"{v:.1f}", ha="center", va="bottom", fontsize=6.0)
    axc.set_xticks(x)
    axc.set_xticklabels(labs)
    axc.set_ylim(0, 145)
    axc.set_yticks([0, 25, 50, 75, 100])
    axc.set_ylabel("Percent (95% Wilson CI)")
    axc.legend(loc="upper right", ncol=1, handlelength=1.0)
    axc.set_title("23-chromosome LOCO\n(n = 300)", fontsize=6.8, pad=4)
    save(fig, "Fig3_comparison_and_LOCO", outdir, formats)


# ==========================================================================
# FIGURE 4  --  structured failure modes
# ==========================================================================
def fig4(outdir, formats):
    fig = plt.figure(figsize=(W2, 72 * MM))
    axa = fig.add_axes([0.285, 0.17, 0.33, 0.72])
    axb = fig.add_axes([0.725, 0.17, 0.26, 0.72])
    fig_label(fig, "a", 0.004, 0.985)
    fig_label(fig, "b", 0.645, 0.985)

    # ---- a: class counts, split by SpliceAI status -----------------------
    cls = FN_CLASSES[::-1]
    y = np.arange(len(cls))
    det = [c[3] for c in cls]
    und = [c[2] - c[3] for c in cls]
    axa.barh(y, det, color=C["blue"], height=0.62, label="SpliceAI \u2265 0.20", lw=0)
    axa.barh(y, und, left=det, color=C["orange"], height=0.62, label="SpliceAI < 0.20", lw=0)
    for yi, c in zip(y, cls):
        axa.text(c[2] + 0.4, yi, f"{c[2]} ({100 * c[2] / 38:.1f}%)", va="center", fontsize=6.0)
    axa.set_yticks(y)
    axa.set_yticklabels([f"{c[0]}  {c[1]}" for c in cls], fontsize=6.0)
    axa.set_xlim(0, 19.5)
    axa.set_xlabel("False negatives (n; total 38)")
    axa.legend(loc="lower right", handlelength=1.0)
    axa.set_title("Mechanistic classes", fontsize=6.8, pad=4, loc="left")
    axa.text(18.9, 6.55, "SpliceAI detected 31/38", fontsize=6.0, ha="right", va="center", color=C["text2"])

    # ---- b: WT score vs |dPWM| ---------------------------------------------
    cmap = {1: C["blue"], 2: C["sky"], 3: C["green"], 4: C["orange"],
            5: C["pink"], 6: C["red"], 7: C["dark"]}
    path = OPTIONAL_DATA["fn_variants"]
    if path.exists():
        rows = read_csv(path)
        for r in rows:
            c = int(float(r["cls"]))
            det_ = int(float(r.get("spliceai_ge_020", 0))) == 1
            axb.scatter(float(r["wt_score"]), abs(float(r["dpwm_splice"])), s=16,
                        facecolor=cmap[c] if det_ else "white", edgecolor=cmap[c], lw=0.9, zorder=3)
        sub = f"{len(rows)} variants (open = SpliceAI < 0.20)"
    else:
        note(f"{path.name} not found; Fig 4b drawn from class medians (reported values). "
             "Add figdata/fn_variants.csv for the per-variant scatter.")
        for c, lab, n, d, wt, dp in FN_CLASSES:
            x0 = 5.0 if wt is None else wt
            axb.scatter(x0, dp, s=22 + 14 * n, facecolor=cmap[c], edgecolor="white", lw=0.5,
                        zorder=3)
            if n >= 5 and c != 2:
                axb.text(x0, dp, str(c), ha="center", va="center", fontsize=6.0,
                         color="white", zorder=4)
            else:                       # small / overlapping markers: label beside, markers stay at true position
                dx, dy = {2: (-0.05, 0.42), 6: (0.42, 0.0)}.get(c, (0.30, 0.0))
                axb.text(x0 + dx, dp + dy, str(c), ha="center", va="center", fontsize=6.0,
                         color=C["dark"], zorder=4)
        axb.annotate("", xy=(4.2, 1.10), xytext=(5.0, 1.10),
                     arrowprops=dict(arrowstyle="-|>", lw=0.6, color=C["orange"]))
        sub = "class medians; marker area scales with n"
    axb.axvline(5.6, ls="--", lw=0.6, color=C["grey"])
    axb.axhline(2.8, ls="--", lw=0.6, color=C["grey"])
    axb.text(5.7, 5.55, "WT \u2265 5.6 bits\n(active donor)", fontsize=6.0, color=C["text2"], va="top")
    axb.text(-0.15, 2.95, "\u03c4_splice = 2.8 bits", fontsize=6.0, color=C["text2"], va="bottom")
    axb.text(0.1, 5.5, "weak wild-type,\nlarge drop", fontsize=6.0, color=C["green"], va="top")
    axb.text(7.4, 0.2, "sub-threshold\nshift", fontsize=6.0, color=C["blue"], ha="right", va="bottom")
    axb.set_xlim(-0.3, 7.6)
    axb.set_ylim(-0.3, 5.7)
    axb.set_xlabel("Wild-type motif score (bits)")
    axb.set_ylabel("|\u0394PWM_splice| (bits)")
    axb.set_title("Wild-type score vs motif shift", fontsize=6.8, pad=4, loc="left")
    axb.text(0.0, -0.27, sub, transform=axb.transAxes, fontsize=6.0, color=C["text2"], ha="left", va="top")
    save(fig, "Fig4_false_negative_classes", outdir, formats)


# ==========================================================================
# FIGURE 5  --  MPRA audit (29 libraries) and TERT promoter
# ==========================================================================
def fig5(outdir, formats):
    p = np.array([l[7] for l in LIBS])
    q = bh_qvalues(p)
    qd = {l[0]: qq for l, qq in zip(LIBS, q)}
    if abs(qd["MSMB"] - 0.025) > 0.002:
        warn(f"MSMB BH q = {qd['MSMB']:.4f}, report says 0.025")
    tq = [qd[k] for k in ("TERT-GAa", "TERT-GBM", "TERT-GSc", "TERT-HEK")]
    if abs(min(tq) - 0.36) > 0.03 or abs(max(tq) - 1.0) > 0.01:
        warn(f"TERT BH q range {min(tq):.2f}-{max(tq):.2f}, report says 0.36-1.00")

    tert = sorted([l for l in LIBS if l[1] == "tert"], key=lambda l: -l[4])
    ind = sorted([l for l in LIBS if l[1] == "ind"], key=lambda l: -l[4])
    rows = [("hdr", "Pooled estimates")] + [("pool", r) for r in POOLED]
    rows += [("gap", "")] + [("hdr", "TERT promoter libraries (in registry)")] + [("lib", r) for r in tert]
    rows += [("gap", "")] + [("hdr", "Independent libraries (25)")] + [("lib", r) for r in ind]

    fig = plt.figure(figsize=(W2, 150 * MM))
    axa = fig.add_axes([0.17, 0.115, 0.20, 0.845])
    axb = fig.add_axes([0.745, 0.625, 0.245, 0.30])
    axc = fig.add_axes([0.805, 0.145, 0.145, 0.265])
    fig_label(fig, "a", 0.004, 0.985)
    fig_label(fig, "b", 0.665, 0.985)
    fig_label(fig, "c", 0.665, 0.50)
    NX = -0.80                       # x (axes units) of the name column

    # ---- a: forest -----------------------------------------------------------
    ytr = axa.get_yaxis_transform()
    for i, (kind, r) in enumerate(rows):
        if kind == "hdr":
            axa.text(NX, i, r, transform=ytr, fontsize=6.4, fontweight="bold",
                     va="center", ha="left", clip_on=False)
            continue
        if kind == "gap":
            continue
        if kind == "pool":
            name, tp, calls, lr, lo, hi, pv = r
            axa.plot([lo, hi], [i, i], color=C["dark"], lw=1.2)
            axa.plot(lr, i, "D", color=C["dark"], ms=4.2)
            qtxt = f"P = {pv:.1e}"
            ctext = f"{tp}/{calls:,}"
        else:
            name, role, tp, calls, lr, lo, hi, pv, ha = r
            col = C["orange"] if role == "tert" else C["blue"]
            sig = qd[name] < 0.05
            axa.plot([lo, hi], [i, i], color=col, lw=0.9)
            axa.plot(lr, i, "s", ms=3.8, mfc=col if sig else "white", mec=col, mew=0.8)
            name = name + ("\u2020" if ha else "")
            qtxt = f"{pv:.3f} / {qd[r[0]]:.2f}" if pv >= 0.001 else f"{pv:.1e} / {qd[r[0]]:.3f}"
            ctext = f"{tp}/{calls}"
        axa.text(NX, i, name, transform=ytr, fontsize=6.2, va="center", ha="left",
                 clip_on=False, color=C["dark"])
        axa.text(1.04, i, f"{lr:.2f} ({lo:.2f}\u2013{hi:.2f})", transform=ytr, fontsize=6.0,
                 va="center", ha="left", clip_on=False)
        axa.text(1.66, i, ctext, transform=ytr, fontsize=6.0, va="center", ha="left", clip_on=False)
        axa.text(2.03, i, qtxt, transform=ytr, fontsize=6.0, va="center", ha="left", clip_on=False)
    axa.text(1.04, -2.2, "LR+ (95% CI)", transform=ytr, fontsize=6, fontweight="bold", ha="left", va="center", clip_on=False)
    axa.text(1.66, -2.2, "Calls TP/n", transform=ytr, fontsize=6, fontweight="bold", ha="left", va="center", clip_on=False)
    axa.text(2.03, -2.2, "P / BH q", transform=ytr, fontsize=6, fontweight="bold", ha="left", va="center", clip_on=False)
    axa.axvline(1, color=C["dark"], lw=0.6, ls="--")
    axa.set_xscale("log")
    axa.set_xlim(0.02, 500)
    axa.set_ylim(len(rows) - 0.4, -3.6)
    axa.set_yticks([])
    axa.spines["left"].set_visible(False)
    axa.set_xticks([0.1, 1, 10, 100])
    axa.set_xticklabels(["0.1", "1", "10", "100"])
    axa.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    axa.set_xlabel("Positive likelihood ratio (LR+, log scale)")
    axa.text(0.85, -2.2, "\u2190 calls depleted for\nfunctional alleles", fontsize=6.0,
             ha="right", va="center", color=C["text2"], clip_on=False, linespacing=1.2)
    axa.text(1.18, -2.2, "calls enriched for\nfunctional alleles \u2192", fontsize=6.0,
             ha="left", va="center", color=C["text2"], clip_on=False, linespacing=1.2)
    fig.text(0.004, 0.012, "Filled = BH q < 0.05;  \u2020 Haldane\u2013Anscombe corrected",
             fontsize=6.0, ha="left", va="bottom", color=C["text2"])

    # ---- b: TERT core promoter (SNV track or summary) ------------------------
    path = OPTIONAL_DATA["tert_snvs"]
    if path.exists():
        data = read_csv(path)
        pos = np.array([float(r["pos"]) for r in data])
        eff = np.array([float(r["mpra_effect"]) for r in data])
        act = np.array([int(float(r["activating"])) for r in data]) == 1
        call = np.array([int(float(r["call_denovo"])) for r in data]) == 1
        axb.scatter(pos[~act], eff[~act], s=4, color=C["mid"], lw=0, label="Neutral")
        axb.scatter(pos[act], eff[act], s=24, color=C["orange"], marker="D", lw=0, label="MPRA-activating", zorder=3)
        axb.scatter(pos[call], eff[call], s=20, facecolor="none", edgecolor=C["red"], lw=0.7,
                    label="NC-MutScan call (de novo mode)", zorder=4)
        axb.axhline(0, color=C["dark"], lw=0.5)
        for xx, yy in zip(pos[act], eff[act]):
            if int(xx) == 1295113:
                axb.annotate("\u2212124 C>T", (xx, yy), xytext=(6, 4), textcoords="offset points", fontsize=6.0)
        axb.set_xlabel("GRCh38 chr5 position")
        axb.set_ylabel("MPRA effect (log2)")
        axb.ticklabel_format(axis="x", style="plain", useOffset=False)
        axb.tick_params(axis="x", labelrotation=30)
        axb.legend(loc="upper left", fontsize=6.0, handlelength=1.0)
        axb.set_title("TERT core promoter, 753 SNVs", fontsize=6.8, pad=4, loc="left")
    else:
        note(f"{path.name} not found; Fig 5b drawn as a summary of reported counts. "
             "Add figdata/tert_snvs.csv for the per-position track.")
        groups = [(m[0], m[1], m[2], m[3], m[4]) for m in TERT_MODES]
        x = np.arange(len(groups))
        w = 0.36
        for j, (metric, col) in enumerate([("Sensitivity (4 activating)", C["blue"]),
                                           ("Specificity (749 neutral)", C["orange"])]):
            vals, lo, hi = [], [], []
            for _, tp, fp, tn, fn in groups:
                k, n = (tp, tp + fn) if j == 0 else (tn, tn + fp)
                lo_, hi_ = wilson(k, n)
                vals.append(100 * k / n)
                lo.append(100 * k / n - 100 * lo_)
                hi.append(100 * hi_ - 100 * k / n)
            xs = x + (j - 0.5) * w
            bar_with_ci(axb, xs, vals, lo, hi, w, col, metric)
            for xx, v, h in zip(xs, vals, hi):
                axb.text(xx, v + h + 1.5, f"{v:.1f}", ha="center", va="bottom", fontsize=6.0)
        axb.set_xticks(x)
        axb.set_xticklabels([g[0] for g in groups], fontsize=6.0)
        axb.set_ylim(0, 135)
        axb.set_yticks([0, 25, 50, 75, 100])
        axb.set_ylabel("Percent (95% Wilson CI)")
        axb.legend(loc="upper left", fontsize=6.0, handlelength=1.0)
        axb.set_title("TERT core promoter, 753 SNVs", fontsize=6.8, pad=4, loc="left")
        axb.text(0.5, -0.30, "SpliceAI: 0/4 (splicing only;\nspecificity not applicable)",
                 transform=axb.transAxes, fontsize=6.0, ha="center", va="top", color=C["text2"])

    # ---- c: pooled MPRA metrics ------------------------------------------------
    metrics = [("Specificity", 95.25, 17427, 18297), ("PPV", 15.29, 157, 1027),
               ("Base rate (tiered)", 18.85, 4250, 22547), ("Sensitivity", 3.69, 157, 4250)]
    ys = np.arange(len(metrics))[::-1]
    for yi, (lab, v, k, n) in zip(ys, metrics):
        lo, hi = wilson(k, n)
        axc.barh(yi, v, color=C["blue"] if lab != "Base rate (tiered)" else C["mid"], height=0.6,
                 lw=0, zorder=2)
        axc.errorbar(v, yi, xerr=[[v - 100 * lo], [100 * hi - v]], fmt="none", ecolor=C["dark"],
                     elinewidth=0.6, capsize=1.5, capthick=0.6, zorder=3)
        axc.text(max(v, 100 * hi) + 2, yi, f"{v:.2f}%", va="center", fontsize=6.0)
    axc.set_yticks(ys)
    axc.set_yticklabels([m[0] for m in metrics])
    axc.set_xlim(0, 112)
    axc.set_xlabel("Pooled over 29 libraries (%)")
    axc.text(-0.95, 1.04, "Pooled MPRA metrics\n(22,547 tiered obs.)", transform=axc.transAxes, fontsize=6.8,
             ha="left", va="bottom", linespacing=1.2)
    axc.text(-0.95, -0.36, "Pooled LR+ = 0.78 (0.66\u20130.92):\ncalls slightly depleted for functional alleles",
             transform=axc.transAxes, fontsize=6.0, ha="left", va="top", color=C["text2"])
    save(fig, "Fig5_MPRA_audit_and_TERT", outdir, formats)


# ==========================================================================
FIGS = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Generate manuscript figures 1-5.")
    ap.add_argument("--only", type=int, nargs="*", choices=sorted(FIGS),
                    help="figure numbers to build (default: all)")
    ap.add_argument("--outdir", type=Path, default=HERE / "paper_figures")
    ap.add_argument("--formats", nargs="*", default=["png", "pdf"],
                    choices=["png", "pdf", "svg", "tiff"])
    ap.add_argument("--dpi", type=int, default=300, help="raster resolution (default 300)")
    args = ap.parse_args(argv)
    set_style()
    mpl.rcParams["savefig.dpi"] = args.dpi
    for n in (args.only or sorted(FIGS)):
        print(f"Figure {n}")
        FIGS[n](args.outdir, args.formats)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())