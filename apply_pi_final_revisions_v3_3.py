"""
apply_pi_final_revisions_v3_3.py
Auto-locates NC-DriverEngine_v3.3_PI_Revised_Manuscript_2.docx (or falls back to
NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx), applies PI Audit Points 1-7,
and outputs NC-DriverEngine_v3.3_PI_Revised_Manuscript_3.docx.
"""

import os
import re
import shutil
import docx
from docx.shared import Pt, RGBColor


def set_para_text(p, text: str, bold=None, size=None):
    fn = "Calibri"
    fs = size or Pt(10.5)
    bd = bold
    if p.runs:
        r0 = p.runs[0]
        if r0.font.name:
            fn = r0.font.name
        if size is None and r0.font.size:
            fs = r0.font.size
        if bd is None:
            bd = bool(r0.bold)
    for r in list(p.runs):
        p._p.remove(r._r)
    nr = p.add_run(text)
    nr.font.name = fn
    nr.font.size = fs
    if bd is not None:
        nr.bold = bd


def ensure_pre_cleanups(doc):
    """Applies the 8 text cleanups in case fallback _1.docx is used."""
    for i, p in enumerate(doc.paragraphs):
        t = p.text
        if not t.strip():
            continue
        if t.startswith("2.9 Leakage-free Leave-One-Out Cross-Validation") and "Random" in t:
            set_para_text(p, "2.9 Leakage-free Leave-One-Out Cross-Validation (LOOCV) classification")
        elif t.startswith("3.7 Divergent driver topologies in melanoma and breast cancer,"):
            set_para_text(p, "3.7 Divergent driver topologies in melanoma and breast cancer")
        elif t.startswith("Table 10. Representative recurrent regulatory loci") and "," in t:
            set_para_text(p, "Table 10. Representative recurrent regulatory loci in melanoma and breast cancer (GRCh38/hg38)")
        elif "M = 56,912,480 non-overlapping 50-bp windows in the independent PCAWG" in t:
            set_para_text(
                p,
                t.replace(
                    "across all M = 56,912,480 non-overlapping 50-bp windows in the independent PCAWG Skin-Melanoma (MELA-AU) WGS cohort (N = 70 donors; 7,987,562 SNVs; donor-level mean μ0 = 0.1342, dispersion r = 0.7289; Table 6).",
                    "across M = 56,912,480 non-N callable 50-bp windows (μ0 = 0.1342, dispersion r = 0.7289; out of M = 61,877,699 total chr1–22,X,Y windows with μ0 = 0.12909, r = 0.6838 in Table 6) in the independent PCAWG Skin-Melanoma (MELA-AU) WGS cohort (N = 70 donors; 7,987,562 SNVs)."
                )
            )
        elif i > 110 and t.startswith("An apparent extreme window at hg19 chr7:140,453,135"):
            set_para_text(
                p,
                "Figure 13. Focal promoter discovery, cis-eQTL comparison, and tissue-specific regulatory topologies. "
                "(A) Local −log10 qlocal across promoter and intronic windows in TCGA-SKCM exomes (prior to excluding the coding BRAF p.V600E exon-15 hotspot at hg19 chr7:140,453,135). "
                "(B) Unadjusted pan-cohort eQTL estimates prior to the tumor-histology and t(14;18) SV adjustments reported in Table 9. "
                "(C) Absolute distance to nearest TSS separating melanoma core-promoter ETS/GABPA hotspots (≤300 bp) from breast cancer intronic enhancers (+15.5 kb to +135.4 kb; Table 10)."
            )
        elif t.startswith("Figure 12.") and "Boxes show median" in t:
            set_para_text(
                p,
                "Figure 12. Expression of BCL6, BCL2, PAX5, and WDR74 in wild-type and promoter-mutated PCAWG donors, "
                "log2(FPKM-UQ + 1). Plotted inset statistics show initial pan-cohort regressions prior to histology fixed effects; "
                "as detailed in Table 9, adding tumor histology fixed effects and adjusting for co-occurring BCL2/t(14;18) "
                "structural variants attenuates BCL6 (qBH = 0.210) and PAX5 (qBH = 0.758) to non-significance while BCL2 retains "
                "a positive expression association within B-cell lymphomas (β = +0.667, p = 2.09 × 10−4)."
            )
        elif "dev_chr1_11.vcf; n = 150, 75 positives and 75 negatives" in t:
            new_t = t.replace(
                "into a development cohort (chr1–chr11; dev_chr1_11.vcf; n = 150, 75 positives and 75 negatives) and a locked blind cohort (chr12–chr22; locked_blind_chr12_22.vcf; n = 150, 75 positives and 75 negatives).",
                "into a calibration pool of N = 955 variants on chr1–11 and chrX (471 Pathogenic/Likely pathogenic and 484 Benign/Likely benign, from which a balanced 150-variant subset dev_chr1_11.vcf was also sampled for v3.1 vs. v3.3 exploratory comparison) and a locked blind test cohort of N = 150 variants on chr12–22 (locked_blind_chr12_22.vcf; 75 Pathogenic/Likely pathogenic and 75 Benign/Likely benign)."
            ).replace(
                "Exactly 75 pathogenic/likely pathogenic and 75 benign/likely benign variants were sampled uniformly without replacement (random.seed(42)) for each chromosome partition (chr1–chr11 and chr12–chr22)",
                "Thresholds were calibrated across the full chr1–11, X calibration set (N = 955), while 75 pathogenic/likely pathogenic and 75 benign/likely benign variants were sampled uniformly without replacement (random.seed(42)) for the locked chr12–22 blind test cohort (N = 150)"
            )
            set_para_text(p, new_t)
        elif "On the development cohort, v3.1 single-sequence scoring reached a sensitivity of 72.0%" in t:
            new_t = t.replace(
                "On the development cohort, v3.1 single-sequence scoring reached a sensitivity of 72.0% but a specificity of only 21.3% (precision 47.8%; F1 score 0.575; Figure 22B). The differential WT-versus-MUT scoring in v3.3 gave up some sensitivity (66.67%) for a large gain in specificity (89.33%), but because its thresholds were calibrated on the development cohort, these development values are in-sample and only the blind values are held out.",
                "On the balanced 150-variant chr1–11 exploratory subset (Figure 22B), v3.1 single-sequence scoring reached 72.0% sensitivity at 21.3% specificity (F1 = 0.575), whereas v3.3 differential scoring reached 66.67% sensitivity and 89.33% specificity in-sample (and 52.87% sensitivity, 95.66% specificity across the full N = 955 chr1–11, X calibration pool; Table 15)."
            )
            set_para_text(p, new_t)
        elif t.startswith("NC-DriverEngine v3.1 takes raw somatic variants through BMR-calibrated hotspot discovery"):
            set_para_text(p, "")

    if len(doc.tables) > 9:
        for row in doc.tables[9].rows[1:]:
            if "2,510,095" in row.cells[1].text:
                row.cells[1].text = "chr16:2,460,094 (hg19:2,510,095)"


def main():
    candidates = [
        "NC-DriverEngine_v3.3_PI_Revised_Manuscript_2.docx",
        os.path.join(os.path.expanduser("~"), "Downloads", "NC-DriverEngine_v3.3_PI_Revised_Manuscript_2.docx"),
        os.path.join(os.path.expanduser("~"), "Desktop", "NC-DriverEngine_v3.3_PI_Revised_Manuscript_2.docx"),
        "NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx",
    ]
    src = next((p for p in candidates if os.path.exists(p)), None)
    if src is None:
        raise FileNotFoundError("Could not find NC-DriverEngine_v3.3_PI_Revised_Manuscript_2.docx or _Manuscript.docx")

    print("Using source file:", src)
    dst = "NC-DriverEngine_v3.3_PI_Revised_Manuscript_3.docx"
    doc = docx.Document(src)
    ensure_pre_cleanups(doc)

    # Identify elements belonging to the 25-locus LOOCV block (Sec 2.9, Sec 3.8, Table 11, Figs 14-15)
    sec29_p = None
    sec29_body_p = None
    sec38_paras = []
    fig14_img_p = None
    fig14_cap_p = None
    fig15_img_p = None
    fig15_cap_p = None
    tbl11 = next((t for t in doc.tables if "LOOCV ROC-AUC" in t.rows[0].cells[2].text), doc.tables[10])

    for idx, p in enumerate(doc.paragraphs):
        t = p.text.strip()
        if t.startswith("2.9 Leakage-free Leave-One-Out Cross-Validation"):
            sec29_p = p
        elif idx > 0 and doc.paragraphs[idx - 1].text.strip().startswith("2.9 Leakage-free Leave-One-Out"):
            sec29_body_p = p
        elif t.startswith("3.8 Leakage-free Leave-One-Out Cross-Validation"):
            sec38_paras.append(p)
        elif t.startswith("Under strict Leave-One-Out Cross-Validation (LOOCV)"):
            sec38_paras.append(p)
        elif t.startswith("Table 11. Leakage-free Leave-One-Out"):
            sec38_paras.append(p)
        elif t.startswith("Table 11 reports exact Average Precision"):
            sec38_paras.append(p)
        elif t.startswith("Figure 14. Out-of-fold Receiver Operating Characteristic"):
            fig14_cap_p = p
            if idx > 0 and not doc.paragraphs[idx - 1].text.strip():
                fig14_img_p = doc.paragraphs[idx - 1]
        elif t.startswith("Figure 15. Out-of-fold ROC and PR curves"):
            fig15_cap_p = p
            if idx > 0 and not doc.paragraphs[idx - 1].text.strip():
                fig15_img_p = doc.paragraphs[idx - 1]

    # Text updates across main body paragraphs
    for i, p in enumerate(doc.paragraphs):
        t = p.text
        if not t.strip():
            continue
        new_t = t

        # Point 2 & 3: Section 1 (3D Hi-C wording + unified v3.3 version naming)
        if new_t.startswith("We built NC-DriverEngine to address both."):
            new_t = (
                "We built NC-DriverEngine v3.3 to address both. It combines statistical hotspot discovery "
                "(BMR-aware sliding windows and zero-included Negative Binomial whole-genome calibration) with a "
                "biophysical reading of each variant (paired wild-type versus mutant ΔPWM log-odds shifts in "
                "NC-MutScan v3.3), adds ENCODE v4 epigenomic context, and filters AID/SHM passengers explicitly "
                "(with 3D Hi-C chromatin-loop integration planned as a future extension). Building on earlier "
                "internal prototypes (v3.0/v3.1), NC-DriverEngine v3.3 integrates the statistical and epigenomic "
                "backbone, the 15-class motif aligner, the variant-overlapping differential ΔPWM scorer, and the "
                "dual-tier API service. We tested the system on 353,368 deduplicated somatic mutations from 466 MAF "
                "donors (470 clinical cases) in the TCGA skin cutaneous melanoma cohort (TCGA-SKCM), 608 PCAWG "
                "whole-genome donors, and 1,105 strictly noncoding ClinVar variants."
            )

        # Point 2 & 3: Section 2.6 (Hi-C as planned + unified v3.3)
        elif new_t.startswith("We mapped candidate loci to the ENCODE cCRE registry"):
            new_t = (
                "We mapped candidate loci to the ENCODE cCRE registry [6] (Release V4) and classified them as "
                "promoter-like (PLS), proximal enhancer-like (pELS), distal enhancer-like (dELS), CTCF insulator "
                "or DNase-only. In NC-DriverEngine v3.3, the lookup runs through the registry's REST API. Linking "
                "mutations in distal elements to target genes via Northwestern 3D Genome Hi-C/Micro-C loop-anchor "
                "data [42–44] is planned for future structural-variant extensions and is not evaluated in the "
                "results below."
            )

        # Point 3 & 5: Renumber Sections 2.10 -> 2.9, 2.11 -> 2.10, 2.12 -> 2.11 and unify v3.3
        elif new_t.startswith("2.10 Dual-mode architecture"):
            new_t = "2.9 Dual-mode architecture and multi-threaded server"
        elif new_t.startswith("NC-DriverEngine v3.1 runs in two synchronized tiers."):
            new_t = (
                new_t.replace("NC-DriverEngine v3.1", "NC-DriverEngine v3.3")
                     .replace("NC-MutScan v3.1", "NC-MutScan v3.3")
                     .replace("21_verify_and_test_pipeline.py", "verify_pipeline.py")
            )
        elif new_t.startswith("Figure 1. Dual-mode architecture"):
            new_t = (
                "Figure 1. Dual-mode architecture of NC-DriverEngine v3.3. Left: cohort-level discovery pipeline "
                "integrating 96-trinucleotide-context BMR modeling, 50-bp sliding-window Poisson scanning, 100-kb "
                "regional BMR calibration, and ENCODE v4 cCRE intersection (3D Hi-C loop mapping is a planned future "
                "module; the 25-locus LOOCV sanity check is detailed in Supplementary Note S1). Right: the real-time "
                "NC-MutScan v3.3 module, which performs sequence sanitization, dual-strand Dirichlet-smoothed ΔPWM "
                "scoring, multi-TF competitive resolution and WRCY/RGYW somatic hypermutation filtering."
            )
        elif new_t.startswith("2.11 NC-MutScan v3.1 biophysical scoring module"):
            new_t = "2.10 NC-MutScan v3.3 biophysical scoring module"
        elif "(Section 3.11)" in new_t or "(Section 3.10)" in new_t or "(Section 3.12)" in new_t or "(Section 2.12)" in new_t:
            new_t = (
                new_t.replace("Section 3.10", "Section 3.9")
                     .replace("Section 3.11", "Section 3.10")
                     .replace("Section 3.12", "Section 3.11")
                     .replace("Section 2.12", "Section 2.11")
                     .replace("Section 2.11", "Section 2.10")
                     .replace("nc_mutscan_heldout_benchmark.py", "benchmark_50_loci.py")
            )
        elif new_t.startswith("2.12 Verification suite and evaluation set"):
            new_t = "2.11 Verification suite and evaluation set"

        # Point 5: Renumber Sections 3.9 -> 3.8, 3.10 -> 3.9, 3.11 -> 3.10, 3.12 -> 3.11 & Tables 12-15 -> 11-14
        elif new_t.startswith("3.9 API performance"):
            new_t = "3.8 API performance"
        elif new_t.startswith("Table 12 and Figure 16"):
            new_t = "Table 11 and Figure 14 report latency and memory for the four service endpoints under concurrent load."
        elif new_t.startswith("Table 12. Multi-threaded API benchmarks"):
            new_t = "Table 11. Multi-threaded API benchmarks"
        elif new_t.startswith("3.10 NC-MutScan v3.1: resolving complex sequence artifacts"):
            new_t = "3.9 NC-MutScan v3.3 (v3.1 baseline suite): resolving complex sequence artifacts"
        elif new_t.startswith("During development we hardened the biophysical module"):
            new_t = new_t.replace("Table 13; Figures 17 and 18", "Table 12; Figures 15 and 16")
        elif "Section 2.11" in new_t:
            new_t = new_t.replace("Section 2.11", "Section 2.10")
        elif new_t.startswith("Table 13. NC-MutScan"):
            new_t = new_t.replace("Table 13.", "Table 12.")
        elif new_t.startswith("3.11 Comparison with conservation-based predictors"):
            new_t = "3.10 Comparison with conservation-based predictors and evaluation-set performance"
        elif new_t.startswith("To assess diagnostic accuracy beyond unit testing"):
            new_t = new_t.replace("Section 3.12", "Section 3.11")
        elif new_t.startswith("To check generalization on real genomic sequence"):
            new_t = (
                new_t.replace("Section 2.12", "Section 2.11")
                     .replace("Table 14; Figure 19", "Table 13; Figure 17")
            )
        elif new_t.startswith("Table 14. 54-locus"):
            new_t = new_t.replace("Table 14.", "Table 13.")
        elif new_t.startswith("On the separate PCAWG core noncoding driver evaluation set"):
            new_t = new_t.replace("Figures 20 and 21", "Figures 18 and 19")
        elif new_t.startswith("3.12 Chromosome-held-out ClinVar benchmark"):
            new_t = "3.11 Chromosome-held-out ClinVar benchmark of strictly noncoding variants"
        elif new_t.startswith("To evaluate NC-MutScan v3.3 on an external clinical benchmark"):
            new_t = new_t.replace("Table 15", "Table 14")
        elif new_t.startswith("Component ablations on the locked chr12–22 set"):
            new_t = new_t.replace("Table 15", "Table 14")
        elif new_t.startswith("Table 15. NC-MutScan v3.3 performance"):
            new_t = new_t.replace("Table 15.", "Table 14.")
        elif new_t.startswith("On the balanced 150-variant chr1–11 exploratory subset"):
            new_t = new_t.replace("Figure 22B", "Figure 20B").replace("Table 15", "Table 14")

        # Point 1: Fix ClinVar comparator contradiction in Section 3.11 (formerly 3.12)
        elif new_t.startswith("The ClinVar cohorts and the PCAWG evaluation set used for CADD"):
            new_t = (
                "CADD v1.6, GWAVA, and FunSeq2 were evaluated only on the separate PCAWG core driver set "
                "(Section 3.10) rather than on the ClinVar cohorts, so Table 14 supports no direct ranking "
                "against those three tools; by contrast, SpliceAI (DS_max ≥ 0.20) was benchmarked directly on "
                "the same locked chr12–22 ClinVar test set (N = 139 scored variants; Table 14)."
            )

        # Renumber Figures 16..22 -> Figures 14..20
        elif new_t.startswith("Figure 16."):
            new_t = new_t.replace("Figure 16.", "Figure 14.").replace("Table 12", "Table 11")
        elif new_t.startswith("Figure 17."):
            new_t = new_t.replace("Figure 17.", "Figure 15.").replace("Section 2.11", "Section 2.10")
        elif new_t.startswith("Figure 18."):
            new_t = new_t.replace("Figure 18.", "Figure 16.")
        elif new_t.startswith("Figure 19."):
            new_t = new_t.replace("Figure 19.", "Figure 17.")
        elif new_t.startswith("Figure 20."):
            new_t = new_t.replace("Figure 20.", "Figure 18.").replace("Table 14", "Table 13")
        elif new_t.startswith("Figure 21."):
            new_t = new_t.replace("Figure 21.", "Figure 19.")
        elif new_t.startswith("Figure 22."):
            new_t = new_t.replace("Figure 22.", "Figure 20.")

        # Point 2, 3, 5 & missing period: Section 4 (Discussion)
        elif new_t.startswith("NC-DriverEngine v3.1 pairs statistical discovery"):
            new_t = (
                "NC-DriverEngine v3.3 pairs statistical discovery with mechanistic interpretation. The BMR-aware "
                "sliding windows separate focal promoter hotspots from UV-driven regional hypermutation, the ENCODE "
                "v4 layer supplies cis-regulatory context, and the ΔPWM module asks why a hotspot matters. In "
                "melanoma, single C>T transitions in core promoters create de novo ETS/GABPA consensus sites [14, 15], "
                "and UV damage is especially efficient at ETS-like sequences [68, 69]. We saw this pattern at the "
                "canonical TERT promoter [12, 13] and again at AC040160.1, AL358472.4/MIR6737 and TCF3. In 433 PCAWG "
                "donors with matched whole-genome and RNA-seq data, unadjusted regression associated mutations at "
                "BCL6, BCL2 and PAX5 with higher expression, whereas the recurrent WDR74 hotspot showed no association. "
                "Incorporating tumor histology fixed effects (+ histology_abbreviation) and within-lymphoma "
                "stratification resolved lineage confounding, showing that apparent expression gains at BCL6 "
                "(qBH = 0.210) and PAX5 (qBH = 0.758) attenuate to non-significance while BCL2 retains a "
                "translocation-co-occurring expression association (β = +0.667, p = 2.09 × 10−4). The two cancers "
                "differed in where their candidates sat: melanoma showed promoter-centric ETS changes and breast "
                "cancer showed distal-enhancer GATA3 loss [66]. The 15-class aligner applies the same logic to "
                "enhancer insertions (TAL1/MYB, LMO2), pioneer-factor elements (FOXA1/AR), CTCF insulators and "
                "noncoding splice sites such as MET exon 14 [8, 36, 40, 41]."
            )
        elif new_t.startswith("Version 3.1 goes after two recurring sources of false calls."):
            new_t = new_t.replace("Version 3.1 goes after", "NC-MutScan v3.3 addresses")
        elif new_t.startswith("Third, the 15-case NC-MutScan verification suite"):
            new_t = (
                new_t.replace("Section 2.11", "Section 2.10")
                     .replace("on the 25-locus benchmark", "on the 25-locus benchmark (Supplementary Note S1)")
            )
        elif new_t.startswith("Fourth, the ClinVar benchmark has its own limits."):
            new_t = (
                new_t.replace("Figure 22A", "Figure 20A")
                     .replace("Section 2.11", "Section 2.10")
            )
        elif new_t.startswith("Several limitations define the scope"):
            new_t = new_t.replace("(Section 3.8)", "(Supplementary Note S1; Ref. [50])")

        # Point 3 & 4: Section 6 (Declarations and data availability)
        elif new_t.startswith("Source code & reproducibility scripts:"):
            new_t = (
                "Source code & reproducibility scripts: https://github.com/mbm25027-stack/noncoding-mutation-platform- "
                "(server.py, sequence_aligner_engine.py, engine_script.js, nc_mutscan_v3_3.py, "
                "build_clinvar_stratified_splits.py, evaluate_clinvar_blind.py, run_pancancer_nb_calibration.py, "
                "run_covariate_eqtl.py, run_loocv_sanity_check.py, verify_pipeline.py, benchmark_50_loci.py, and "
                "generate_publication_figures.py)."
            )
        elif new_t.startswith("Data: TCGA-SKCM somatic mutation data"):
            new_t = (
                "Data: TCGA-SKCM somatic mutation data are publicly available from the NCI Genomic Data Commons "
                "(Release 46.0) [28, 29]. PCAWG consensus somatic variants, RNA-seq expression (FPKM-UQ), segmented "
                "copy number, consensus purity, and histology tables are publicly available from the UCSC Xena PCAWG "
                "hub [49, 59]. Locked statistical verification outputs are stored in results/tables/ "
                "(covariate_adjusted_eqtl_results.csv, pancancer_nb_bmr_summary.csv, "
                "clinvar_chr12_22_blind_results.csv, and benchmark_hotspot_matrix_25loci.csv), and 300-DPI "
                "publication figures are stored in results/figures/. Noncoding ClinVar variants [72] were drawn "
                "from the GRCh38 ClinVar VCF (accessed November 2026); the calibration set (dev_chr1_11.vcf), "
                "locked blind cohort (locked_blind_chr12_22.vcf), and benchmark scripts are archived in "
                "https://github.com/mbm25027-stack/noncoding-mutation-platform-."
            )
        elif new_t.startswith("Author contributions:"):
            new_t = (
                "Author contributions: S.B. (lead researcher) and A.B. (principal investigator) conceived the study. "
                "S.B. designed the NC-DriverEngine v3.3 framework (including the NC-MutScan v3.3 biophysical scoring "
                "module), performed the bioinformatics benchmarking and statistical analyses, and wrote the manuscript. "
                "A.B. supervised the project and revised the manuscript. Both authors approved the final manuscript."
            )

        if new_t != t:
            set_para_text(p, new_t)

    # Point 5: Relocate the 25-locus LOOCV section, Table 11, and Figures 14-15 to Supplementary Note S1 at end
    body = doc.element.body
    supp_hdr = doc.add_paragraph("Supplementary Note S1. Internal LOOCV calibration sanity check (25 curated loci)")
    supp_hdr.runs[0].bold = True
    supp_hdr.runs[0].font.name = "Calibri"
    supp_hdr.runs[0].font.size = Pt(12)
    supp_hdr.runs[0].font.color.rgb = RGBColor(0x1E, 0x3A, 0x5F)

    if sec29_p is not None:
        set_para_text(sec29_p, "S1.1 Leakage-free Leave-One-Out Cross-Validation (LOOCV) methodology", bold=True, size=Pt(10.5))
        body.append(sec29_p._p)
    if sec29_body_p is not None:
        set_para_text(
            sec29_body_p,
            sec29_body_p.text.replace("38_repaired_ml_hotspot_matrix.csv", "benchmark_hotspot_matrix_25loci.csv")
        )
        body.append(sec29_body_p._p)

    for sp in sec38_paras:
        st = sp.text
        if st.startswith("3.8 "):
            set_para_text(sp, "S1.2 Out-of-fold LOOCV sanity-check results", bold=True, size=Pt(10.5))
        elif st.startswith("Under strict Leave-One-Out"):
            set_para_text(
                sp,
                st.replace("38_repaired_ml_hotspot_matrix.csv", "benchmark_hotspot_matrix_25loci.csv")
                  .replace("Table 11, Figures 14 and 15", "Supplementary Table S1 and Supplementary Figures S1–S2")
            )
        elif st.startswith("Table 11. "):
            set_para_text(sp, st.replace("Table 11.", "Supplementary Table S1."), bold=True, size=Pt(9.5))
        elif st.startswith("Table 11 reports"):
            set_para_text(
                sp,
                st.replace("Table 11", "Supplementary Table S1")
                  .replace("Figures 14–15", "Supplementary Figures S1–S2"),
                size=Pt(9)
            )
        body.append(sp._p)

    for idx, p in enumerate(doc.paragraphs):
        if p.text.startswith("Supplementary Table S1."):
            p._p.addnext(tbl11._tbl)
            break

    if fig14_img_p is not None:
        body.append(fig14_img_p._p)
    if fig14_cap_p is not None:
        set_para_text(fig14_cap_p, fig14_cap_p.text.replace("Figure 14.", "Supplementary Figure S1."), size=Pt(9))
        body.append(fig14_cap_p._p)

    if fig15_img_p is not None:
        body.append(fig15_img_p._p)
    if fig15_cap_p is not None:
        set_para_text(fig15_cap_p, fig15_cap_p.text.replace("Figure 15.", "Supplementary Figure S2."), size=Pt(9))
        body.append(fig15_cap_p._p)

    doc.save(dst)
    print("Saved fully revised manuscript:", dst)

    # Cleanly copy/alias internal script filenames in workspace so repo matches manuscript
    renames = {
        "21_verify_and_test_pipeline.py": "verify_pipeline.py",
        "nc_mutscan_heldout_benchmark.py": "benchmark_50_loci.py",
        "22_nb_bmr_and_covariate_eqtl.py": "run_covariate_eqtl.py",
        "23_day3_pancancer_nb_and_ml_benchmark.py": "run_pancancer_nb_calibration.py",
        "25_day4_repaired_ml_loocv_and_figures.py": "run_loocv_sanity_check.py",
    }
    for old_f, new_f in renames.items():
        if os.path.exists(old_f) and not os.path.exists(new_f):
            shutil.copy2(old_f, new_f)
            print(f"Created clean script alias: {old_f} -> {new_f}")


if __name__ == "__main__":
    main()
