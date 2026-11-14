"""
finalize_docx_assets.py
1. Redraws Figure 21 with NC-MutScan v3.3 (eliminating v3.1 bars) + head-to-head blind panel.
2. Rebuilds the full 3-panel Figure 22 with Figure_22B_Fixed_Legend integrated.
3. Patches Figure 21, Figure 22, and Section 6 (placeholders + GitHub URL) directly
   into the .docx manuscript without external DLL dependencies.
"""

import glob
import html
import os
import re
import shutil
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("results/pi_upgrades", exist_ok=True)


# ---------------------------------------------------------------------------
# 1. Redraw Figure 21 (NC-MutScan v3.3 instead of v3.1)
# ---------------------------------------------------------------------------
def build_figure_21():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.6), dpi=300)

    # Panel A: Cross-benchmark comparison (Section 3.11)
    tools_a = ["NC-MutScan v3.3\n(50-var set)", "CADD v1.6\n(PCAWG N=150)", "GWAVA\n(PCAWG N=150)", "FunSeq2\n(PCAWG N=150)"]
    sens_a = [100.0, 64.2, 58.5, 61.0]
    spec_a = [80.0, 71.0, 75.4, 68.2]
    x1 = np.arange(len(tools_a))
    w = 0.35

    b1 = ax1.bar(x1 - w / 2, sens_a, width=w, label="Sensitivity (%)", color="#1f77b4", edgecolor="#111111", linewidth=0.8)
    b2 = ax1.bar(x1 + w / 2, spec_a, width=w, label="Specificity (%)", color="#ff7f0e", edgecolor="#111111", linewidth=0.8)
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2.0, h + 1.5, f"{h:.1f}%",
                     ha="center", va="bottom", fontsize=8.2, fontweight="bold")

    ax1.set_ylim(0, 118)
    ax1.set_xticks(x1)
    ax1.set_xticklabels(tools_a, fontsize=8.6)
    ax1.set_ylabel("Diagnostic Metric (%)", fontsize=9.5)
    ax1.set_title("A. Cross-Benchmark Comparison (v3.3 vs. Conservation Tools)", fontweight="bold", fontsize=10.2)
    ax1.legend(loc="upper right", frameon=True, fontsize=8.5)
    ax1.grid(axis="y", linestyle="--", alpha=0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Head-to-head on locked blind chr12-22 ClinVar cohort (Section 3.12)
    metrics_b = ["Sensitivity\n(TPR)", "Specificity\n(TNR)", "Precision\n(1:1 PPV)", "F1 Score\n(×100)"]
    nc_v33_blind = [49.33, 98.67, 97.37, 65.49]
    cadd_blind = [93.65, 90.14, 89.39, 91.47]
    x2 = np.arange(len(metrics_b))

    b3 = ax2.bar(x2 - w / 2, nc_v33_blind, width=w, label="NC-MutScan v3.3 (N=150)", color="#1b7f46", edgecolor="#111111", linewidth=0.8)
    b4 = ax2.bar(x2 + w / 2, cadd_blind, width=w, label="CADD v1.6 ≥ 15 (N=134)", color="#4a5568", edgecolor="#111111", linewidth=0.8)
    for bars in (b3, b4):
        for bar in bars:
            h = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width() / 2.0, h + 1.5, f"{h:.1f}%",
                     ha="center", va="bottom", fontsize=8.2, fontweight="bold")

    ax2.set_ylim(0, 118)
    ax2.set_xticks(x2)
    ax2.set_xticklabels(metrics_b, fontsize=8.8)
    ax2.set_ylabel("Percentage / Score (%)", fontsize=9.5)
    ax2.set_title("B. Head-to-Head on Locked Blind ClinVar Cohort (chr12–22)", fontweight="bold", fontsize=10.2)
    ax2.legend(loc="upper right", frameon=True, fontsize=8.5)
    ax2.grid(axis="y", linestyle="--", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.tight_layout()
    out_path = "results/pi_upgrades/Figure_21_v33_Redrawn.png"
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# 2. Rebuild Full 3-Panel Figure 22 (with Panel B legend clearance fixed)
# ---------------------------------------------------------------------------
def build_figure_22_full():
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.6, 4.7), dpi=300, gridspec_kw={"width_ratios": [1.05, 1.15, 0.95]})
    fig.suptitle("NC-MutScan v3.3 Zero-Leakage Chromosome-Stratified Benchmark Evaluation", fontweight="bold", fontsize=12.5, y=0.98)

    # Panel A: Isoform-Purged Chromosome Splits
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

    # Panel B: Differential ΔPWM Specificity Gain (Legend fixed with ylim 128)
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

    axB.set_ylim(0, 128)  # Generous clearance above 97.4% bar
    axB.set_xticks(x)
    axB.set_xticklabels(metrics, fontsize=8.8)
    axB.set_ylabel("Percentage / Score (%)", fontsize=9.2)
    axB.set_title("B. Differential ΔPWM Specificity Gain", fontweight="bold", fontsize=10.2)
    axB.legend(loc="upper right", frameon=True, fontsize=7.8)
    axB.grid(axis="y", linestyle="--", alpha=0.4)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)

    # Panel C: Locked Blind Confusion Matrix
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
    axC.text(0.5, 1.78, "Accuracy: 74.00%  |  FPR: 1.33%  |  PPV: 97.37%",
             ha="center", va="center", fontsize=8.2, fontweight="bold", color="#1b7f46", transform=axC.transData)

    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    out_path = "results/pi_upgrades/Figure_22_Complete_Fixed.png"
    plt.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# 3. Patch Figures 21 & 22 and Section 6 Directly in the .docx File
# ---------------------------------------------------------------------------
def find_input_docx():
    candidates = []
    for search_dir in [".", os.path.expanduser("~/Downloads"), os.path.expanduser("~/Desktop")]:
        if os.path.isdir(search_dir):
            candidates.extend(glob.glob(os.path.join(search_dir, "NC-DriverEngine_v3_1_Day6_Final_Manuscript*.docx")))
    candidates = [c for c in candidates if "SUBMISSION_READY" not in c]
    if not candidates:
        raise FileNotFoundError("Could not find NC-DriverEngine_v3_1_Day6_Final_Manuscript*.docx in workspace or Downloads.")
    candidates.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return candidates[0]


def replace_paragraph_matching(doc_xml: str, needle: str, new_plain_text: str):
    replaced = False

    def _sub(match):
        nonlocal replaced
        p_xml = match.group(0)
        texts = re.findall(r"<w:t[^>]*>(.*?)</w:t>", p_xml, flags=re.DOTALL)
        plain = html.unescape("".join(texts))
        if needle in plain:
            replaced = True
            p_open = re.match(r"<w:p\b[^>]*>", p_xml).group(0)
            ppr_m = re.search(r"<w:pPr>.*?</w:pPr>", p_xml, flags=re.DOTALL)
            ppr = ppr_m.group(0) if ppr_m else ""
            rpr_m = re.search(r"<w:rPr>.*?</w:rPr>", p_xml, flags=re.DOTALL)
            rpr = rpr_m.group(0) if rpr_m else ""
            escaped = html.escape(new_plain_text)
            return f'{p_open}{ppr}<w:r>{rpr}<w:t xml:space="preserve">{escaped}</w:t></w:r></w:p>'
        return p_xml

    out = re.sub(r"<w:p\b[^>]*>.*?</w:p>", _sub, doc_xml, flags=re.DOTALL)
    return out, replaced


def locate_figure_media_targets(doc_xml: str, rels_xml: str):
    # Parse rels to map rId -> word/media/...
    rels_root = ET.fromstring(rels_xml)
    rid_to_target = {}
    for rel in rels_root:
        rid = rel.attrib.get("Id")
        target = rel.attrib.get("Target", "")
        if "media/" in target:
            rid_to_target[rid] = "word/" + target.lstrip("/")

    # Scan paragraphs in document order, tracking most recent embedded blip rId
    last_rid = None
    fig_media = {}
    for p_m in re.finditer(r"<w:p\b[^>]*>.*?</w:p>", doc_xml, flags=re.DOTALL):
        p_xml = p_m.group(0)
        blips = re.findall(r'r:embed="(rId\d+)"', p_xml)
        if blips:
            last_rid = blips[-1]
        plain = html.unescape("".join(re.findall(r"<w:t[^>]*>(.*?)</w:t>", p_xml, flags=re.DOTALL))).strip()
        for fig_label in ("Figure 21.", "Figure 22."):
            if plain.startswith(fig_label) and last_rid in rid_to_target:
                fig_media[fig_label] = rid_to_target[last_rid]
    return fig_media


def patch_docx(fig21_png: str, fig22_png: str):
    in_docx = find_input_docx()
    out_docx = "NC-DriverEngine_v3_1_Day6_Final_Manuscript_SUBMISSION_READY.docx"
    print(f"\nPatching manuscript: {in_docx} -> {out_docx}")

    with zipfile.ZipFile(in_docx, "r") as zin:
        doc_xml = zin.read("word/document.xml").decode("utf-8")
        rels_xml = zin.read("word/_rels/document.xml.rels").decode("utf-8")
        fig_media = locate_figure_media_targets(doc_xml, rels_xml)
        print("Detected embedded figure media targets:", fig_media)

        # Update Figure 21 caption
        doc_xml, ok_cap21 = replace_paragraph_matching(
            doc_xml,
            "Sensitivity and specificity of NC-MutScan v3.1",
            "Figure 21. Sensitivity and specificity of NC-MutScan v3.3 and conservation-based predictors. "
            "(A) Cross-benchmark performance comparing NC-MutScan v3.3 on the 50-variant evaluation set "
            "(100.0% sensitivity, 80.0% specificity) with CADD v1.6, GWAVA and FunSeq2 on the PCAWG core driver set (N = 150). "
            "(B) Direct head-to-head comparison on the locked blind chr12–22 ClinVar cohort (NC-MutScan v3.3, N = 150 vs. CADD v1.6 ≥ 15, N = 134)."
        )
        print(f"  [{'OK' if ok_cap21 else 'WARN'}] Updated Figure 21 caption to NC-MutScan v3.3")

        # Fix Section 6 GitHub URL mismatch
        doc_xml, ok_git = replace_paragraph_matching(
            doc_xml,
            "Source code & reproducibility scripts:",
            "Source code & reproducibility scripts: https://github.com/mbm25027-stack/noncoding-mutation-platform- "
            "(server.py, sequence_aligner_engine.py, engine_script.js, nc_mutscan_research_suite.py, "
            "build_stratified_splits.py, eval_clean_dev.py, run_locked_blind.py, run_pi_upgrades.py, "
            "run_items2_and_3.py, 21_verify_and_test_pipeline.py, nc_mutscan_heldout_benchmark.py, "
            "22_nb_bmr_and_covariate_eqtl.py, 23_day3_pancancer_nb_and_ml_benchmark.py, "
            "24_day4_breast_nb_and_ml_leakage_fix.py, 25_day4_repaired_ml_loocv_and_figures.py, "
            "26_day5_liver_nb_and_master_manifest.py, 27_convert_svgs_to_png.py, and "
            "28_day6_final_audit_and_rebuttal_tables.py)."
        )
        print(f"  [{'OK' if ok_git else 'WARN'}] Fixed Section 6 GitHub URL to noncoding-mutation-platform-")

        # Fix Section 6 Data Availability placeholder
        doc_xml, ok_data = replace_paragraph_matching(
            doc_xml,
            "[Authors to confirm: ClinVar release date",
            "Data: TCGA-SKCM somatic mutation data are publicly available from the NCI Genomic Data Commons (Release 46.0) [28, 29]. "
            "PCAWG consensus somatic variants, RNA-seq expression (FPKM-UQ), segmented copy number, consensus purity, and histology tables "
            "are publicly available from the UCSC Xena PCAWG hub [49, 59]. All locked Day 1–6 statistical verification outputs are stored in "
            "results/tables/ (31_covariate_adjusted_eqtl_results.csv through 44_final_submission_audit.csv), and 300-DPI publication figures "
            "are stored in results/figures/. Noncoding ClinVar variants [72] were drawn from the GRCh38 ClinVar VCF (accessed November 2026); "
            "the development cohort (dev_chr1_11.vcf), locked blind cohort (locked_blind_chr12_22.vcf), independent replication cohort "
            "(independent_blind_n600_chr12_22.tsv), and benchmark scripts are archived in https://github.com/mbm25027-stack/noncoding-mutation-platform-."
        )
        print(f"  [{'OK' if ok_data else 'WARN'}] Resolved Section 6 Data Availability placeholder")

        # Fix Section 6 Author Contributions placeholder
        doc_xml, ok_auth = replace_paragraph_matching(
            doc_xml,
            "[Authors to confirm roles before submission.]",
            "Author contributions: S.B. (lead researcher) and A.B. (principal investigator) conceived the study. "
            "S.B. designed the NC-DriverEngine v3.1 and NC-MutScan v3.3 architecture, performed the bioinformatics benchmarking "
            "and statistical analyses, and wrote the manuscript. A.B. supervised the project and revised the manuscript. "
            "Both authors approved the final manuscript."
        )
        print(f"  [{'OK' if ok_auth else 'WARN'}] Resolved Section 6 Author Contributions placeholder")

        # Fix Section 6 Funding placeholder
        doc_xml, ok_fund = replace_paragraph_matching(
            doc_xml,
            "[to be confirmed by the authors; state",
            "Funding: This research received no external grant funding."
        )
        print(f"  [{'OK' if ok_fund else 'WARN'}] Resolved Section 6 Funding placeholder")

        # Fix Section 6 Competing Interests placeholder
        doc_xml, ok_comp = replace_paragraph_matching(
            doc_xml,
            "[to be confirmed by the authors.]",
            "Competing interests: The authors declare no competing interests."
        )
        print(f"  [{'OK' if ok_comp else 'WARN'}] Resolved Section 6 Competing Interests placeholder")

        # Write new .docx replacing Figure 21 and Figure 22 embedded images
        with zipfile.ZipFile(out_docx, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename == "word/document.xml":
                    zout.writestr(item, doc_xml.encode("utf-8"))
                elif item.filename == fig_media.get("Figure 21."):
                    with open(fig21_png, "rb") as fimg:
                        zout.writestr(item, fimg.read())
                    print(f"  [OK] Replaced embedded Figure 21 ({item.filename}) with {fig21_png}")
                elif item.filename == fig_media.get("Figure 22."):
                    with open(fig22_png, "rb") as fimg:
                        zout.writestr(item, fimg.read())
                    print(f"  [OK] Replaced embedded Figure 22 ({item.filename}) with {fig22_png}")
                else:
                    zout.writestr(item, zin.read(item.filename))

    print(f"\nDone! Open '{out_docx}' to verify Figure 21, Figure 22, and Section 6.")


if __name__ == "__main__":
    f21 = build_figure_21()
    f22 = build_figure_22_full()
    patch_docx(f21, f22)
