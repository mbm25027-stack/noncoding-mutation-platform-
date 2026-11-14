import re
import docx
from docx.shared import Pt

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

path = "NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx"
doc = docx.Document(path)

for i, p in enumerate(doc.paragraphs):
    t = p.text
    if not t.strip():
        continue

    # 1. Clean appended clauses from Section 2.9, Section 3.7 headings & Table 10 caption
    if t.startswith("2.9 Leakage-free Leave-One-Out Cross-Validation"):
        set_para_text(p, "2.9 Leakage-free Leave-One-Out Cross-Validation (LOOCV) classification")
    elif t.startswith("3.7 Divergent driver topologies in melanoma and breast cancer"):
        set_para_text(p, "3.7 Divergent driver topologies in melanoma and breast cancer")
    elif t.startswith("Table 10. Representative recurrent regulatory loci"):
        set_para_text(p, "Table 10. Representative recurrent regulatory loci in melanoma and breast cancer (GRCh38/hg38)")

    # 2. Reconcile callable vs total 50-bp windows in Section 3.3b
    elif "M = 56,912,480 non-overlapping 50-bp windows" in t:
        set_para_text(
            p,
            t.replace(
                "across all M = 56,912,480 non-overlapping 50-bp windows in the independent PCAWG Skin-Melanoma (MELA-AU) WGS cohort (N = 70 donors; 7,987,562 SNVs; donor-level mean μ0 = 0.1342, dispersion r = 0.7289; Table 6).",
                "across M = 56,912,480 non-N callable 50-bp windows (μ0 = 0.1342, dispersion r = 0.7289; out of M = 61,877,699 total chr1–22,X,Y windows with μ0 = 0.12909, r = 0.6838 in Table 6) in the independent PCAWG Skin-Melanoma (MELA-AU) WGS cohort (N = 70 donors; 7,987,562 SNVs)."
            )
        )

    # 3. Restore Figure 13 caption where duplicate BRAF paragraph overwrote it (after Table 10)
    elif i > 110 and t.startswith("An apparent extreme window at hg19 chr7:140,453,135"):
        set_para_text(
            p,
            "Figure 13. Focal promoter discovery, cis-eQTL comparison, and tissue-specific regulatory topologies. "
            "(A) Local −log10 qlocal across promoter and intronic windows in TCGA-SKCM exomes (prior to excluding the coding BRAF p.V600E exon-15 hotspot at hg19 chr7:140,453,135). "
            "(B) Unadjusted pan-cohort eQTL estimates prior to the tumor-histology and t(14;18) SV adjustments reported in Table 9. "
            "(C) Absolute distance to nearest TSS separating melanoma core-promoter ETS/GABPA hotspots (≤300 bp) from breast cancer intronic enhancers (+15.5 kb to +135.4 kb; Table 10)."
        )

    # 4. Clarify Figure 12 caption regarding unadjusted panels vs Table 9 adjusted models
    elif t.startswith("Figure 12."):
        set_para_text(
            p,
            "Figure 12. Expression of BCL6, BCL2, PAX5, and WDR74 in wild-type and promoter-mutated PCAWG donors, "
            "log2(FPKM-UQ + 1). Plotted inset statistics show initial pan-cohort regressions prior to histology fixed effects; "
            "as detailed in Table 9, adding tumor histology fixed effects and adjusting for co-occurring BCL2/t(14;18) "
            "structural variants attenuates BCL6 (qBH = 0.210) and PAX5 (qBH = 0.758) to non-significance while BCL2 retains "
            "a positive expression association within B-cell lymphomas (β = +0.667, p = 2.09 × 10−4)."
        )

    # 5. Reconcile Section 2.12 & Section 3.12 ClinVar N=955 vs n=150 descriptions
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

    # 6. Fix the broken sentence fragment in Section 4 (Discussion)
    elif ".66 × 10−3 to 2.47 × 10−5)" in t:
        new_t = re.sub(
            r"\.66\s*×\s*10−3\s*to\s*2\.47\s*×\s*10−5\)\s*and\s*PAX5\s*\(3\.86\s*×\s*10−2\s*to\s*4\.29\s*×\s*10−3\)\s*associations,\s*consistent\s*with\s*removal\s*of\s*variance\s*from\s*copy-number\s*changes\s*and\s*stromal\s*dilution\.",
            "",
            t
        )
        set_para_text(p, new_t)

    # 7. Remove duplicate old v3.1 Conclusion paragraph
    elif t.startswith("NC-DriverEngine v3.1 takes raw somatic variants through BMR-calibrated hotspot discovery"):
        set_para_text(p, "")

    # 8. Fix residual string artifacts in Section 6 (Declarations)
    elif "pipeline build–6" in t or "blind benchmark (N = 150)_chr12_22.tsv" in t:
        new_t = (
            t.replace("All locked pipeline build–6 statistical verification outputs", "All locked statistical verification outputs")
             .replace(", independent replication cohort (locked chromosome-held-out blind benchmark (N = 150)_chr12_22.tsv),", ",")
        )
        set_para_text(p, new_t)

# Fix C16orf59 GRCh38/hg19 coordinate in Table 10 (table index 9)
for row in doc.tables[9].rows[1:]:
    if "2,510,095" in row.cells[1].text:
        row.cells[1].text = "chr16:2,460,094 (hg19:2,510,095)"

doc.save(path)
print("Successfully cleaned all 8 items in", path)
