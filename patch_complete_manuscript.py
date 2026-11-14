"""
patch_complete_manuscript.py
Synchronizes all text, tables (Table 9 95% CIs, Table 15B CADD head-to-head,
Table 16 38-FN mechanisms), Section 3.13 (N=600 replication cohort),
Abstract/Discussion/Conclusion, Reference [70], Figure 13 (3-panel with 95% CIs),
and Figure 22 (clean Panel C footer spacing) into the final .docx.
"""

import html
import os
import re
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("results/pi_upgrades", exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Rebuild Complete 3-Panel Figure 13 (with 95% CI Forest Plot in Panel B)
# ---------------------------------------------------------------------------
def build_figure_13_complete():
    fig, (axA, axB, axC) = plt.subplots(
        3, 1, figsize=(8.4, 10.6), dpi=300, gridspec_kw={"height_ratios": [1.0, 0.95, 1.15]}
    )

    # Panel A: Local BMR calibration
    loci_a = ["TERT", "AC040160.1", "AL358472.4", "TCF3", "MKRN1", "VRK3", "AC093525.2", "OXNAD1"]
    vals_a = [30.113, 23.851, 19.148, 18.008, 60.641, 36.633, 33.006, 33.006]
    cols_a = ["#b84a27"] * 4 + ["#6b7280"] * 4
    xa = np.arange(len(loci_a))
    axA.bar(xa, vals_a, color=cols_a, width=0.62)
    axA.set_xticks(xa)
    axA.set_xticklabels(loci_a, rotation=28, ha="right", fontsize=8.8)
    axA.set_ylabel(r"$-\log_{10}\, q_{\mathrm{local}}$", fontsize=9.5)
    axA.set_ylim(0, 68)
    axA.set_title("A  Local BMR calibration: focal promoter drivers vs. regional intronic peaks",
                  loc="left", fontweight="bold", fontsize=10.5)
    from matplotlib.patches import Patch
    axA.legend(
        handles=[
            Patch(facecolor="#b84a27", label="Promoter (≤300 bp from TSS)"),
            Patch(facecolor="#6b7280", label="Intronic regional peak"),
        ],
        loc="upper right", frameon=False, fontsize=8.5
    )
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)

    # Panel B: True Forest Plot with 95% Confidence Intervals
    genes_b = ["BCL6", "BCL2", "PAX5", "WDR74\n(neg. control)"]
    betas = [0.903, 1.322, 3.197, -0.158]
    ci_low = [0.653, 0.720, 1.011, -0.690]
    ci_high = [1.153, 1.924, 5.383, 0.374]
    q_labels = [
        "β = +0.903 [+0.653, +1.153]; q = 1.68 × 10⁻¹¹",
        "β = +1.322 [+0.720, +1.924]; q = 4.95 × 10⁻⁵",
        "β = +3.197 [+1.011, +5.383]; q = 5.72 × 10⁻³",
        "β = −0.158 [−0.690, +0.374]; q = 0.559",
    ]
    cols_b = ["#b84a27", "#b84a27", "#b84a27", "#6b7280"]
    yb = np.arange(len(genes_b))[::-1]
    axB.axvline(0, color="#222222", linewidth=1.0, zorder=1)
    for y, b, lo, hi, lbl, col in zip(yb, betas, ci_low, ci_high, q_labels, cols_b):
        axB.errorbar(
            b, y, xerr=[[b - lo], [hi - b]],
            fmt="o", color=col, ecolor=col, elinewidth=2.2, capsize=4, markersize=7, zorder=3
        )
        axB.text(hi + 0.14, y, lbl, va="center", fontsize=8.4, color="#111111")
    axB.set_yticks(yb)
    axB.set_yticklabels(genes_b, fontsize=9.0, fontstyle="italic")
    axB.set_xlim(-1.0, 9.2)
    axB.set_xlabel(r"Covariate-adjusted $\beta_{\mathrm{mut}}$ ($\log_2$ FPKM-UQ + 1, 95% CI)", fontsize=9.5)
    axB.set_title("B  Covariate-adjusted eQTL effects with 95% CIs (matched PCAWG WGS + RNA-seq)",
                  loc="left", fontweight="bold", fontsize=10.5)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    # Panel C: Divergent Regulatory Topologies
    genes_c = ["RPL13A", "CTA-397H3.3", "TRMT10C", "C16orf59", "GATA3", "AKT1", "SF3B1", "AFF2"]
    dists = [117, 10, 36, 14, 15866, 15537, 32981, 135413]
    cols_c = ["#1e4976"] * 4 + ["#b84a27"] * 4
    yc = np.arange(len(genes_c))[::-1]
    axC.axvline(300, color="#6b7280", linestyle="--", linewidth=1.1)
    axC.text(330, 7.15, "300 bp", color="#6b7280", fontsize=8.2)
    for y, d, col in zip(yc, dists, cols_c):
        axC.hlines(y, 5, d, color=col, linewidth=2.0)
        axC.scatter([d], [y], color=col, s=42, zorder=3)
    axC.set_xscale("log")
    axC.set_xlim(5, 3e5)
    axC.set_yticks(yc)
    axC.set_yticklabels(genes_c, fontsize=8.8)
    axC.set_xlabel("Absolute distance to nearest TSS (bp, log scale)", fontsize=9.5)
    axC.set_title("C  Divergent regulatory topologies", loc="left", fontweight="bold", fontsize=10.5)
    axC.legend(
        handles=[
            Patch(facecolor="#1e4976", label="Melanoma core promoter (ETS/GABPA shifts)"),
            Patch(facecolor="#b84a27", label="Breast intronic enhancer"),
        ],
        loc="center right", frameon=False, fontsize=8.5
    )
    axC.spines["top"].set_visible(False)
    axC.spines["right"].set_visible(False)

    plt.tight_layout()
    out_path = "results/pi_upgrades/Figure_13_Complete_CI.png"
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# 2. Rebuild Figure 22 with Clean Panel C Footer Clearance
# ---------------------------------------------------------------------------
def build_figure_22_clean():
    fig, (axA, axB, axC) = plt.subplots(
        1, 3, figsize=(15.6, 4.9), dpi=300, gridspec_kw={"width_ratios": [1.05, 1.15, 0.95]}
    )
    fig.suptitle("NC-MutScan v3.3 Zero-Leakage Chromosome-Stratified Benchmark Evaluation",
                 fontweight="bold", fontsize=12.5, y=0.97)

    # Panel A
    categories = ["Blind Holdout\n(chr12–22)", "Dev Sample\n(chr1–11)", "Blind Pool\n(chr12–22)", "Dev Pool\n(chr1–11)"]
    path_counts = [75, 75, 3472, 6036]
    ben_counts = [75, 75, 193490, 314319]
    y = np.arange(len(categories))
    h = 0.34
    axA.barh(y + h / 2, path_counts, height=h, color="#b83227", label="Pathogenic (Non-Coding)")
    axA.barh(y - h / 2, ben_counts, height=h, color="#2b6cb0", label="Benign (Non-Coding)")
    axA.set_xscale("log")
    axA.set_xlim(10, 2.5e6)
    axA.set_yticks(y)
    axA.set_yticklabels(categories, fontsize=8.8, fontweight="bold")
    axA.set_xlabel("Variant Count (Log Scale)", fontsize=9.2)
    axA.set_title("A. Isoform-Purged Chromosome Splits", fontweight="bold", fontsize=10.2)
    for yi, pv, bv in zip(y, path_counts, ben_counts):
        axA.text(pv * 1.22, yi + h / 2, f"{pv:,}", va="center", fontsize=7.8, color="#8b1e16", fontweight="bold")
        axA.text(bv * 1.22, yi - h / 2, f"{bv:,}", va="center", fontsize=7.8, color="#1d4ed8", fontweight="bold")
    axA.legend(loc="lower right", frameon=True, fontsize=7.8)
    axA.grid(axis="x", linestyle="--", alpha=0.4)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)

    # Panel B
    metrics = ["Sensitivity\n(TPR)", "Specificity\n(TNR)", "Precision\n(PPV)", "F1-Score\n(×100)"]
    v31_dev = [72.0, 21.3, 47.8, 57.5]
    v33_dev = [66.7, 89.3, 86.2, 75.2]
    v33_blind = [49.3, 98.7, 97.4, 65.5]
    x = np.arange(len(metrics))
    w = 0.26
    b1 = axB.bar(x - w, v31_dev, width=w, label="v3.1 Single-Seq (Dev)", color="#9aa0a6")
    b2 = axB.bar(x, v33_dev, width=w, label="v3.3 ΔPWM (Dev chr1–11)", color="#2b6cb0")
    b3 = axB.bar(x + w, v33_blind, width=w, label="v3.3 ΔPWM (Blind chr12–22)", color="#1b7f46")
    for bars in (b1, b2, b3):
        for bar in bars:
            bh = bar.get_height()
            axB.text(bar.get_x() + bar.get_width() / 2.0, bh + 1.5, f"{bh:.1f}",
                     ha="center", va="bottom", fontsize=7.6, fontweight="bold")
    axB.set_ylim(0, 128)
    axB.set_xticks(x)
    axB.set_xticklabels(metrics, fontsize=8.8)
    axB.set_ylabel("Percentage / Score (%)", fontsize=9.2)
    axB.set_title("B. Differential ΔPWM Specificity Gain", fontweight="bold", fontsize=10.2)
    axB.legend(loc="upper right", frameon=True, fontsize=7.8)
    axB.grid(axis="y", linestyle="--", alpha=0.4)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    # Panel C (with clean spacing below x-axis tick labels)
    cm = np.array([[37, 38], [1, 74]])
    cell_colors = np.array([["#7fb3d5", "#85c1e9"], ["#f4f6f7", "#1f618d"]])
    axC.set_xlim(-0.5, 1.5)
    axC.set_ylim(1.5, -0.5)
    for r in range(2):
        for c in range(2):
            rect = plt.Rectangle((c - 0.5, r - 0.5), 1, 1, facecolor=cell_colors[r, c], edgecolor="#222222", linewidth=1.0)
            axC.add_patch(rect)
            pct = cm[r, c] / 75.0 * 100.0
            txt_col = "white" if (r == 1 and c == 1) else "#111111"
            tag = ["TP", "FN", "FP", "TN"][r * 2 + c]
            axC.text(c, r, f"{tag} = {cm[r, c]}\n({pct:.1f}%)", ha="center", va="center",
                     fontsize=9.5, fontweight="bold", color=txt_col)
    axC.set_xticks([0, 1])
    axC.set_xticklabels(["Predicted\nDriver", "Predicted\nBenign"], fontsize=8.8, fontweight="bold")
    axC.set_yticks([0, 1])
    axC.set_yticklabels(["ClinVar\nPathogenic", "ClinVar\nBenign"], fontsize=8.8, fontweight="bold")
    axC.set_title("C. Locked Blind Matrix (chr12–22)", fontweight="bold", fontsize=10.2)
    axC.text(0.5, -0.22, "Accuracy: 74.00%  |  FPR: 1.33%  |  PPV: 97.37%",
             ha="center", va="top", fontsize=8.1, fontweight="bold", color="#1b7f46", transform=axC.transAxes)

    plt.tight_layout(rect=[0, 0.06, 1, 0.92])
    out_path = "results/pi_upgrades/Figure_22_Complete_Fixed.png"
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# 3. XML Patching Helpers for Word (.docx)
# ---------------------------------------------------------------------------
def make_paragraph_xml(text: str, bold: bool = False, space_after: str = "120") -> str:
    escaped = html.escape(text)
    b_tag = "<w:b/>" if bold else ""
    return (
        f'<w:p><w:pPr><w:spacing w:after="{space_after}" w:line="276" w:lineRule="auto"/></w:pPr>'
        f'<w:r><w:rPr>{b_tag}<w:sz w:val="22"/><w:szCs w:val="22"/></w:rPr>'
        f'<w:t xml:space="preserve">{escaped}</w:t></w:r></w:p>'
    )


def make_table_xml(headers: list, rows: list) -> str:
    def _cell(val: str, is_header: bool = False) -> str:
        esc = html.escape(str(val))
        shd = '<w:shd w:val="clear" w:color="auto" w:fill="E8EEF5"/>' if is_header else ""
        b = "<w:b/>" if is_header else ""
        return (
            f'<w:tc><w:tcPr>{shd}<w:tcMar><w:top w:w="70" w:type="dxa"/><w:bottom w:w="70" w:type="dxa"/>'
            f'<w:left w:w="90" w:type="dxa"/><w:right w:w="90" w:type="dxa"/></w:tcMar></w:tcPr>'
            f'<w:p><w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr>'
            f'<w:r><w:rPr>{b}<w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr>'
            f'<w:t xml:space="preserve">{esc}</w:t></w:r></w:p></w:tc>'
        )

    tr_head = "<w:tr>" + "".join(_cell(h, True) for h in headers) + "</w:tr>"
    tr_rows = "".join("<w:tr>" + "".join(_cell(c, False) for c in r) + "</w:tr>" for r in rows)
    borders = (
        '<w:tblBorders>'
        '<w:top w:val="single" w:sz="6" w:space="0" w:color="444444"/>'
        '<w:bottom w:val="single" w:sz="6" w:space="0" w:color="444444"/>'
        '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="CCCCCC"/>'
        '<w:left w:val="none"/><w:right w:val="none"/><w:insideV w:val="none"/>'
        '</w:tblBorders>'
    )
    return f'<w:tbl><w:tblPr><w:tblW w:w="0" w:type="auto"/>{borders}</w:tblPr>{tr_head}{tr_rows}</w:tbl>'


def replace_paragraph_xml(doc_xml: str, needle: str, replacement_xml: str):
    replaced = False

    def _sub(match):
        nonlocal replaced
        p_xml = match.group(0)
        plain = html.unescape("".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p_xml, flags=re.DOTALL)))
        if needle in plain:
            replaced = True
            return replacement_xml
        return p_xml

    out = re.sub(r"<w:p\b[^>]*>.*?</w:p>", _sub, doc_xml, flags=re.DOTALL)
    return out, replaced


def locate_figure_media_targets(doc_xml: str, rels_xml: str):
    rels_root = ET.fromstring(rels_xml)
    rid_to_target = {}
    for rel in rels_root:
        rid = rel.attrib.get("Id")
        target = rel.attrib.get("Target", "")
        if "media/" in target:
            rid_to_target[rid] = "word/" + target.lstrip("/")

    last_rid = None
    fig_media = {}
    for p_m in re.finditer(r"<w:p\b[^>]*>.*?</w:p>", doc_xml, flags=re.DOTALL):
        p_xml = p_m.group(0)
        blips = re.findall(r'r:embed="(rId\d+)"', p_xml)
        if blips:
            last_rid = blips[-1]
        plain = html.unescape("".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p_xml, flags=re.DOTALL))).strip()
        for fig_label in ("Figure 13.", "Figure 21.", "Figure 22."):
            if plain.startswith(fig_label) and last_rid in rid_to_target:
                fig_media[fig_label] = rid_to_target[last_rid]
    return fig_media


# ---------------------------------------------------------------------------
# 4. Apply All Text, Table, and Figure Updates to Create Final .docx
# ---------------------------------------------------------------------------
def build_final_manuscript():
    in_docx = "NC-DriverEngine_v3_1_Day6_Final_Manuscript_SUBMISSION_READY.docx"
    if not os.path.exists(in_docx):
        raise FileNotFoundError(f"Expected {in_docx} in current directory.")
    out_docx = "NC-DriverEngine_v3_1_Final_Submission_Manuscript.docx"

    f13_png = build_figure_13_complete()
    f22_png = build_figure_22_clean()

    with zipfile.ZipFile(in_docx, "r") as zin:
        doc_xml = zin.read("word/document.xml").decode("utf-8")
        rels_xml = zin.read("word/_rels/document.xml.rels").decode("utf-8")
        fig_media = locate_figure_media_targets(doc_xml, rels_xml)
        print("Embedded media targets:", fig_media)

        # A. Update Abstract closing sentences to include N=600 replication & FN mechanism
        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "Distinguishing functional noncoding driver mutations from passengers is hardest",
            make_paragraph_xml(
                "Distinguishing functional noncoding driver mutations from passengers is hardest in tumors with extreme background mutation, "
                "such as UV-driven melanoma and AID-driven lymphoma. We present NC-DriverEngine v3.1, a dual-tier framework combining "
                "sliding-window scanning, local background-rate calibration, zero-included negative binomial (NB) modeling and NC-MutScan, "
                "a module that scores transcription factor binding shifts (ΔPWM) and filters AID/SHM passengers. In 467 TCGA-SKCM melanomas "
                "it recovered the TERT promoter driver (qlocal = 7.71 × 10−31) and three further core-promoter candidates. NB calibration "
                "across 608 PCAWG genomes confirmed GATA3 enhancer recurrence in breast cancer (qNB = 2.74 × 10−5). In 433 PCAWG donors with "
                "matched RNA-seq, mutations at BCL6, BCL2 and PAX5 were associated with higher expression after covariate adjustment "
                "(qeQTL ≤ 5.72 × 10−3). On two independent chromosome-held-out ClinVar benchmarks of strictly noncoding variants (chr12–chr22; "
                "N = 150 and N = 600), frozen NC-MutScan v3.3 achieved 98.67% and 96.33% specificity (97.37% and 93.82% balanced precision) "
                "at 49.33% and 55.67% sensitivity, with 94.7% of missed pathogenic variants localized to deep-intronic or sub-threshold "
                "splice-region sites outside local motif windows."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Updated Abstract")

        # B. Cite Reference [70] in Section 3.3b
        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "Genome-wide NB recalibration in breast cancer confirmed the canonical GATA3",
            make_paragraph_xml(
                "In breast adenocarcinoma (BRCA, N = 120, 770,911 SNVs), the mean rate dropped >10-fold to μ0 = 0.01270 SNVs/50 bp "
                "(σ² = 0.01295, Fano factor = 1.02) while maintaining stable dispersion (α̂ = 1.5407, r = 0.6491). Genome-wide NB "
                "recalibration in breast cancer confirmed the canonical GATA3 intronic enhancer hotspot (chr10:8,111,433, 8 mutated donors, "
                "FE = 364.8, qNB = 2.74 × 10−5) [66, 70] as the sole genome-wide significant noncoding driver while properly attenuating "
                "low-recurrence 2–4 donor loci (AKT1, SF3B1, PLEKHS1, and AFF2; all qNB = 1.00)."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Cited Reference [70] in Section 3.3b")

        # C. Update Table 9 (Section 3.6) with 95% Confidence Intervals & Figure 13 caption
        table9_new = make_table_xml(
            ["Gene", "Matched cohorts", "WT / MUT", "Mean WT vs. MUT", "Unadj. p", "CovAdj βmut (95% CI)", "CovAdj p", "qeQTL (BH)"],
            [
                ["BCL6", "Kidney-RCC, Lymph-BNHL, Lymph-CLL, Ovary-AdenoCA, Panc-AdenoCA", "394 / 96", "3.562 vs. 4.623", "3.58 × 10−15", "+0.903 [+0.653, +1.153]", "4.20 × 10−12", "1.68 × 10−11"],
                ["BCL2", "Lymph-BNHL, Lymph-CLL", "134 / 46", "4.373 vs. 5.298", "2.66 × 10−3", "+1.322 [+0.720, +1.924]", "2.47 × 10−5", "4.95 × 10−5"],
                ["PAX5", "Kidney-RCC, Lymph-BNHL, Lymph-CLL", "291 / 25", "2.756 vs. 5.025", "3.86 × 10−2", "+3.197 [+1.011, +5.383]", "4.29 × 10−3", "5.72 × 10−3"],
                ["WDR74 (neg. control)", "Lymph-BNHL, Ovary-AdenoCA, Panc-AdenoCA", "275 / 5", "3.208 vs. 3.100", "0.693", "−0.158 [−0.690, +0.374]", "0.559", "0.559"],
            ]
        )
        # Replace old Table 9 (the 4th <w:tbl> or the one containing BCL6 and WDR74)
        def _replace_tbl9(m):
            t_xml = m.group(0)
            if "BCL6" in t_xml and "WDR74" in t_xml:
                return table9_new
            return t_xml
        doc_xml, n_t9 = re.subn(r"<w:tbl\b[^>]*>.*?</w:tbl>", _replace_tbl9, doc_xml, flags=re.DOTALL)
        print(f"  [{'OK' if n_t9 else 'WARN'}] Updated Table 9 with 95% CIs")

        # D. Clarify +616.6 kb in Section 3.7 & update Figure 13 caption
        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "TCGA-BRCA candidates [65] instead mapped to distal intronic enhancers",
            make_paragraph_xml(
                "To compare regulatory architecture across tumor types without pre-filtering on expression, we examined the top "
                "motif-impacted noncoding hotspots in each cohort. Melanoma candidates clustered in core promoters (≤300 bp from the TSS) "
                "and showed ETS/GABPA motif changes [15]. Several of them (RPL13A, CTA-397H3.3, TRMT10C and C16orf59) recur in 12–17 donors "
                "and match the classic UV-induced pyrimidine-dimer hotspot pattern at ETS motifs [15, 68, 69]; we therefore treat them as "
                "UV motif-recurrence loci rather than expression-supported drivers. TCGA-BRCA candidates [65] instead mapped to distal "
                "intronic enhancers (+15.5 kb to +616.6 kb across the broader candidate set; four representative loci up to +135.4 kb listed "
                "in Table 10), with GATA3 motif loss at the GATA3 locus [66, 70], a gene that is itself recurrently mutated in breast cancer [67]. "
                "Consistent with the WGS Negative Binomial calibration in Section 3.3b, only the 8-donor GATA3 intronic enhancer hotspot "
                "retained genome-wide significance (qNB = 2.74 × 10−5), whereas the 3–4 donor breast loci (AKT1, SF3B1, and AFF2) attenuated "
                "to qNB = 1.00 and are retained in Table 10 strictly to illustrate tissue-specific enhancer distance distributions."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Clarified +616.6 kb in Section 3.7")

        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "Figure 13. Focal promoter driver discovery, transcriptional validation",
            make_paragraph_xml(
                "Figure 13. Focal promoter driver discovery, transcriptional validation and tissue-specific regulatory topologies. "
                "(A) Local 100-kb BMR calibration separates four recurrent promoter drivers (TERT, AC040160.1, AL358472.4 and TCF3; ≤300 bp from the TSS) "
                "from regional intronic peaks (MKRN1, VRK3/RPL13A, AC093525.2/C16orf59, OXNAD1) in TCGA-SKCM. "
                "(B) Covariate-adjusted cis-eQTL regression coefficients (βmut) with 95% confidence intervals in matched PCAWG donors (Table 9). "
                "(C) Distance to the nearest TSS (log scale) for melanoma core-promoter candidates with ETS/GABPA motif shifts and breast intronic "
                "enhancer candidates (Table 10)."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Updated Figure 13 caption")

        # E. Reconcile Section 3.11 with Figure 21A (Unified v3.3 on 50-var and 54-locus sets)
        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "To assess diagnostic accuracy beyond unit testing, we evaluated NC-MutScan v3.1",
            make_paragraph_xml(
                "To assess diagnostic accuracy under a unified software version, we also evaluated NC-MutScan v3.3 in unpaired single-sequence "
                "fallback mode (wt_seq = None) on the curated 50-variant set (25 canonical promoter and enhancer driver sequences spanning TERT, "
                "SDHD, TAL1, PLEKHS1, and PAX5 motifs; 25 GC-rich promoter and neutral intergenic controls). NC-MutScan v3.3 achieved 100.0% "
                "sensitivity (25/25 true positives; 95% Wilson CI 86.7–100.0%), 80.0% specificity (20/25 true negatives; 95% CI 60.9–91.1%), "
                "83.3% precision, and 90.0% overall accuracy (45/50; Figure 21A), where 5 GC-rich promoter controls matched endogenous Sp1/ETS "
                "motifs in the absence of a paired wild-type reference allele."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Reconciled Section 3.11 (50-var set) with Figure 21A")

        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "To check generalization on real genomic sequence, with no matrix changes",
            make_paragraph_xml(
                "Evaluating the 54-locus GRCh38/hg38 stress-test panel (4 canonical driver hotspots plus 50 unselected GRCh38 intergenic windows) "
                "in unpaired single-sequence mode (wt_seq = None) under NC-MutScan v3.3 recovered 4/4 (100.0%) canonical driver hotspots "
                "(TERT −124 C>T, PLEKHS1, SDHD, and TAL1), isolated 4/50 (8.0%) intergenic windows containing reference N bases as Ambiguous "
                "IUPAC Sequence, classified 9/50 (18.0%) as Benign Polymorphism, and called 37/50 (74.0%) unmodified human intergenic windows "
                "as motif-bearing (41/54 [75.9%] total; Table 14). Because unmodified intergenic DNA frequently harbors pre-existing endogenous "
                "TF consensus motifs, unpaired single-sequence scanning cannot distinguish pre-existing binding sites from de novo somatic motif "
                "creation—directly motivating the paired wild-type versus mutant differential ΔPWM mode evaluated in Sections 3.12–3.13."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Reconciled Section 3.11 (54-locus set) with v3.3")

        # Update Table 14 to reflect v3.3 unpaired single-sequence stress test
        table14_new = make_table_xml(
            ["Evaluation category", "Loci (n)", "NC-MutScan v3.3 unpaired call (wt_seq = None)", "Share", "Mechanistic interpretation"],
            [
                ["Canonical driver hotspots", "4", "4/4 Driver Variant (100.0% recall)", "7.4%", "TERT, PLEKHS1, SDHD, TAL1 motifs detected"],
                ["Unselected hg38 intergenic DNA", "37", "Driver Variant (endogenous motif match)", "68.5%", "Pre-existing endogenous TF motifs; motivates paired ΔPWM"],
                ["Unselected hg38 intergenic DNA", "9", "Benign Polymorphism", "16.7%", "No high-affinity endogenous motif match"],
                ["Reference windows with N bases", "4", "Ambiguous IUPAC Sequence", "7.4%", "Safely sanitized without runtime errors"],
                ["Total", "54", "41 Driver / 9 Benign / 4 Sanitized", "100.0%", "Demonstrates necessity of paired WT-vs.-MUT ΔPWM"],
            ]
        )
        def _replace_tbl14(m):
            t_xml = m.group(0)
            if "Canonical driver hotspots" in t_xml and "54" in t_xml:
                return table14_new
            return t_xml
        doc_xml, n_t14 = re.subn(r"<w:tbl\b[^>]*>.*?</w:tbl>", _replace_tbl14, doc_xml, flags=re.DOTALL)
        print(f"  [{'OK' if n_t14 else 'WARN'}] Updated Table 14 for v3.3")

        # F. Update Section 3.12 (Bayes PPV + Table 16 38-FN Breakdown + CADD v1.6 Head-to-Head + Section 3.13 N=600)
        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "On the development cohort, v3.1 single-sequence scoring reached a sensitivity",
            make_paragraph_xml(
                "On the development cohort, v3.1 single-sequence scoring reached a sensitivity of 72.0% but a specificity of only 21.3% "
                "(precision 47.8%; F1 score 0.575; Figure 22B). The differential WT-versus-MUT scoring in v3.3 gave up some sensitivity (66.67%) "
                "for a large gain in specificity (89.33%), though because thresholds were calibrated on chr1–chr11, these development metrics are "
                "in-sample. On the locked blind cohort (chr12–chr22), specificity rose to 98.67% (74/75; 1.33% false-positive rate) at 49.33% "
                "sensitivity (37/75) and 97.37% balanced precision (37/38). Because the benchmark was balanced 1:1 by design, Bayes' theorem "
                "yields a positive predictive value of 39.90% at the exact 1.76% empirical pathogenic prevalence of the source chr12–chr22 "
                "noncoding ClinVar pool (3,472 of 196,962 variants) and 43.02% at a 2.0% clinical prior (a >21-fold enrichment over baseline "
                "prevalence, with 98.96–99.09% negative predictive value)."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Added Bayes realistic-prevalence PPV to Section 3.12")

        table16_xml = (
            make_paragraph_xml(
                "Table 16. Mechanistic breakdown of the 38 false-negative pathogenic variants in the locked blind chr12–chr22 cohort, "
                "annotated by ClinVar molecular consequence.",
                bold=True, space_after="60"
            )
            + make_table_xml(
                ["ClinVar molecular consequence category", "Missed pathogenic variants (n)", "Share of FNs (%)", "Mechanistic boundary"],
                [
                    ["Deep-intronic (>20 bp from exon junction)", "24", "63.2%", "Outside local promoter/TF PWMs and 9/13-mer splice windows"],
                    ["Splice-region / canonical splice (below cutoff)", "12", "31.6%", "Sub-threshold ΔPWM or non-canonical branch-point/ESE perturbation"],
                    ["Other / intergenic regulatory", "2", "5.3%", "Distal regulatory motifs outside the 9 core TF PWMs"],
                    ["Total blind false negatives", "38", "100.0%", "94.7% (36/38) are deep-intronic or sub-threshold splice variants"],
                ]
            )
            + make_paragraph_xml(
                "Annotating all 38 missed pathogenic blind variants against ClinVar molecular-consequence records (Table 16) showed that "
                "24 (63.2%) were deep-intronic variants (>20 bp from exon boundaries), 12 (31.6%) were splice-region or non-canonical splice "
                "variants whose ΔPWM shifts fell below the specificity-calibrated cutoffs, and 2 (5.3%) were intergenic regulatory variants "
                "outside the nine core TF PWMs. Thus, 94.7% (36/38) of false negatives reflect deep-intronic and sub-threshold splicing "
                "mechanisms outside the local motif scope of NC-MutScan v3.3. Direct head-to-head scoring of the exact same locked blind "
                "cohort with precomputed CADD v1.6 (PHRED ≥ 15; 134/150 variants indexed, excluding 16 unindexed multi-base indels; Figure 21B) "
                "yielded higher sensitivity (93.65%, 59/63; 95% CI 84.8–97.5%) on these largely Mendelian disease-gene variants, whereas "
                "NC-MutScan v3.3 reduced the false-positive rate more than seven-fold (1.33% [1/75] vs. 9.86% [7/71]; specificity 98.67% "
                "vs. 90.14%; balanced precision 97.37% vs. 89.39%) while scoring multi-base indels natively and reporting the exact disrupted TF or splice motif."
            )
            + make_paragraph_xml(
                "3.13 Independent replication on a second disjoint held-out cohort (N = 600, chr12–chr22)",
                bold=True, space_after="80"
            )
            + make_paragraph_xml(
                "To obtain tighter confidence intervals without touching the original 150-variant locked blind set, we sampled a second, "
                "completely non-overlapping cohort of 600 strictly noncoding GRCh38 ClinVar variants from chr12–chr22 (300 pathogenic/likely "
                "pathogenic, 300 benign/likely benign; random.seed(104)), excluding all coding consequences and every coordinate present in "
                "either dev_chr1_11.vcf or locked_blind_chr12_22.vcf. Evaluated once with frozen NC-MutScan v3.3, the classifier achieved "
                "55.67% sensitivity (167/300; 95% CI 50.0–61.2%), 96.33% specificity (289/300; 95% CI 93.6–97.9%; false-positive rate 3.67%, "
                "11/300), 93.82% balanced precision (167/178; 95% CI 89.3–96.5%), 76.00% accuracy (95% CI 72.4–79.2%), and F1 = 0.6987. Under "
                "Bayes' theorem, these operating characteristics correspond to a positive predictive value of 21.41% at the exact 1.76% "
                "empirical pathogenic prevalence of the source chr12–chr22 ClinVar pool (99.18% NPV) and 23.65% at a 2.0% clinical prior "
                "(99.07% NPV), confirming reproducible high-specificity variant prioritization across independent chromosome-split cohorts."
            )
        )
        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "The ClinVar cohorts and the PCAWG evaluation set used for CADD, GWAVA and FunSeq2",
            table16_xml
        )
        print(f"  [{'OK' if ok else 'WARN'}] Inserted Table 16, CADD head-to-head, and Section 3.13 (N=600)")

        # G. Update Section 4 (Limitations, Fourth) & Section 5 (Conclusion)
        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "Fourth, the ClinVar benchmark has its own limits. Sensitivity was 49.33%",
            make_paragraph_xml(
                "Fourth, the ClinVar benchmark clarifies the biophysical scope of NC-MutScan v3.3. Across the two independent chr12–chr22 "
                "held-out cohorts (N = 150 and N = 600), sensitivity was 49.33% (95% CI 38.3–60.4%) and 55.67% (50.0–61.2%) at 98.67% and "
                "96.33% specificity, and mechanistic annotation of the 38 false negatives in the primary blind cohort confirmed that 94.7% "
                "(36/38) were deep-intronic (63.2%) or sub-threshold splice-region variants (31.6%) outside the nine core TF PWMs and local "
                "9/13-mer splice windows. Under realistic prior prevalences (1.76–2.0%), Bayes' theorem yields a positive predictive value of "
                "21.4–43.0% (a 12- to 24-fold enrichment over baseline prevalence, with >98.9% NPV), making NC-MutScan v3.3 well suited as a "
                "high-specificity triage filter complementary to genome-wide conservation scores such as CADD v1.6. Finally, ClinVar pathogenic "
                "and benign classifications largely reflect Mendelian disease-gene variants rather than somatic cancer drivers, and evaluating "
                "NC-MutScan v3.3 in unpaired single-sequence mode (wt_seq = None) on the 54-locus GRCh38 set confirmed that paired WT-versus-MUT "
                "differential scoring is essential to suppress false positives from endogenous intergenic TF motifs."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Updated Section 4 (Limitations, Fourth)")

        doc_xml, ok = replace_paragraph_xml(
            doc_xml,
            "NC-DriverEngine v3.1 takes raw somatic variants through BMR-calibrated hotspot discovery",
            make_paragraph_xml(
                "NC-DriverEngine v3.1 takes raw somatic variants through BMR-calibrated hotspot discovery, zero-included Negative Binomial "
                "dispersion modeling across WGS cohorts, epigenomic and 3D annotation, ΔPWM scoring, and AID/SHM filtering, serving results "
                "through a real-time HTTP microservice and a zero-network browser engine. In TCGA-SKCM it recovered the canonical TERT promoter "
                "driver and nominated three further core-promoter candidates. Across 608 PCAWG WGS donors, zero-included NB calibration "
                "maintained stable dispersion across >10-fold shifts in mutation rate, confirming GATA3 enhancer recurrence in breast cancer "
                "and canonical CTNNB1, KRAS, and TP53 hotspots in hepatocellular carcinoma. In 433 matched PCAWG WGS and RNA-seq donors, "
                "mutations at BCL6, BCL2, and PAX5 significantly increased target gene expression after covariate adjustment, whereas WDR74 "
                "showed no effect. Under leakage-free LOOCV, multi-feature prioritization achieved PR-AUC = 0.8894 (pperm = 0.0050) for driver "
                "prioritization and PR-AUC = 0.8986 (pperm = 0.0030) for scale-invariant architecture classification. Across two independent "
                "chromosome-held-out ClinVar benchmarks of strictly noncoding variants (chr12–chr22; N = 150 and N = 600), frozen NC-MutScan v3.3 "
                "reached 98.67% and 96.33% specificity (97.37% and 93.82% balanced precision) at 49.33% and 55.67% sensitivity."
            )
        )
        print(f"  [{'OK' if ok else 'WARN'}] Updated Section 5 (Conclusion)")

        # H. Write final .docx replacing Figure 13 and Figure 22 images
        with zipfile.ZipFile(out_docx, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename == "word/document.xml":
                    zout.writestr(item, doc_xml.encode("utf-8"))
                elif item.filename == fig_media.get("Figure 13."):
                    with open(f13_png, "rb") as fimg:
                        zout.writestr(item, fimg.read())
                    print(f"  [OK] Replaced Figure 13 ({item.filename}) with 3-panel 95% CI version")
                elif item.filename == fig_media.get("Figure 22."):
                    with open(f22_png, "rb") as fimg:
                        zout.writestr(item, fimg.read())
                    print(f"  [OK] Replaced Figure 22 ({item.filename}) with clean Panel C footer")
                else:
                    zout.writestr(item, zin.read(item.filename))

    print(f"\n[SUCCESS] Generated '{out_docx}' — ready for submission!")


if __name__ == "__main__":
    build_final_manuscript()
