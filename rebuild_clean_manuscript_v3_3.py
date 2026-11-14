"""
rebuild_clean_manuscript_v3_3.py
Performs a clean, single-pass surgical rebuild starting from the unpatched source
(NC-DriverEngine_v3_1_Day6_Final_Manuscript_SUBMISSION_READY.docx) so there are:
  - Zero [Revision Note] tags or post-Declaration pasted blocks.
  - Exact preservation of v3.0 / v3.1 / v3.3 module history (Section 2.7 = v3.0,
    Section 2.11 & 3.10 = v3.1 verification/54-locus benchmark, Section 3.12 = v3.3 ClinVar).
  - Exact Ensembl-verified GRCh38 coordinates across all tables and prose, purging
    hg19 chr7:140,453,135 (BRAF p.V600E) from noncoding tables.
  - Reframed contribution: transparent statistical & biophysical triage filter, UV-hotspot
    detector (TERT confirmed driver; AL358472.4, AC040160.1, TCF3 reclassified as recurrent
    UV-associated ETS-motif hotspots), and interpretable delta-PWM tool.
  - Honest PCAWG cis-eQTL reporting (Table 9): BCL6 (q=0.210) and PAX5 (q=0.758) attenuate
    after histology adjustment; BCL2 retains association (beta=+0.667, p=2.09e-4 after
    t(14;18)/SV adjustment) with explicit disclosure that 42/46 mutated donors carry BCL2 SVs.
  - Honest ClinVar & Saturation-Mutagenesis benchmarking (Table 15): SpliceAI (F1=0.933) vs.
    NC-MutScan v3.3 (F1=0.655), ablations showing splice-PWM dominance on ClinVar, and
    TERT 251-bp core-promoter saturation mutagenesis (N=753 SNVs).
  - Resolution of all Review 1 open items: title reframed ("Sequence-to-Regulatory-Topology"
    without claiming unanalyzed Hi-C), near-Poisson NB Fano factors (1.02-1.19) acknowledged,
    25-locus LOOCV framed as internal calibration check, TTN labeled strictly as negative
    control in Table 7, References 50 & 70 cited, and "Repaired"/"Day 1-6"/"n600" scrubbed.
"""

import os
import re
import docx
from docx.shared import Pt, RGBColor
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn


def set_cell_shading(cell, hex_color: str):
    shading = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{hex_color}"/>')
    cell._tc.get_or_add_tcPr().append(shading)


def set_cell_margins(cell, top=55, bottom=55, left=75, right=75):
    tcPr = cell._tc.get_or_add_tcPr()
    tcMar = OxmlElement("w:tcMar")
    for side, val in (("top", top), ("bottom", bottom), ("left", left), ("right", right)):
        node = OxmlElement(f"w:{side}")
        node.set(qn("w:w"), str(val))
        node.set(qn("w:type"), "dxa")
        tcMar.append(node)
    tcPr.append(tcMar)


def set_table_borders(table, border_color="CBD5E1"):
    tblPr = table._tbl.tblPr
    borders = parse_xml(
        f'<w:tblBorders {nsdecls("w")}>'
        f'<w:top w:val="single" w:sz="6" w:space="0" w:color="{border_color}"/>'
        f'<w:bottom w:val="single" w:sz="8" w:space="0" w:color="1E3A5F"/>'
        f'<w:insideH w:val="single" w:sz="4" w:space="0" w:color="{border_color}"/>'
        f'<w:left w:val="none"/><w:right w:val="none"/><w:insideV w:val="none"/>'
        f'</w:tblBorders>'
    )
    tblPr.append(borders)


def replace_table_in_place(doc, old_tbl, headers: list[str], rows: list[list[str]], highlight_keys=()):
    new_tbl = doc.add_table(rows=1 + len(rows), cols=len(headers))
    new_tbl.autofit = True
    set_table_borders(new_tbl)

    trPr = new_tbl.rows[0]._tr.get_or_add_trPr()
    trPr.append(OxmlElement("w:tblHeader"))
    trPr.append(OxmlElement("w:cantSplit"))
    for j, h_text in enumerate(headers):
        c = new_tbl.rows[0].cells[j]
        set_cell_shading(c, "1E3A5F")
        set_cell_margins(c, top=70, bottom=70, left=75, right=75)
        p = c.paragraphs[0]
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(0)
        r = p.add_run(h_text)
        r.bold = True
        r.font.name = "Calibri"
        r.font.size = Pt(8.5)
        r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    for i, rvals in enumerate(rows):
        row = new_tbl.rows[i + 1]
        row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
        is_hi = any(k in rvals[0] for k in highlight_keys)
        bg = "EFF6FF" if is_hi else ("F8FAFC" if i % 2 == 1 else "FFFFFF")
        for j, val in enumerate(rvals):
            c = row.cells[j]
            set_cell_shading(c, bg)
            set_cell_margins(c, top=50, bottom=50, left=75, right=75)
            p = c.paragraphs[0]
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(0)
            r = p.add_run(val)
            r.font.name = "Calibri"
            r.font.size = Pt(8.5)
            r.font.color.rgb = RGBColor(0x0F, 0x17, 0x2A)
            if j == 0:
                r.bold = True

    old_tbl._element.addprevious(new_tbl._element)
    old_tbl._element.getparent().remove(old_tbl._element)
    return new_tbl


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


def main():
    src = "NC-DriverEngine_v3_1_Day6_Final_Manuscript_SUBMISSION_READY.docx"
    dst = "NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx"
    doc = docx.Document(src)

    # -----------------------------------------------------------------------
    # 1. Surgical Paragraph Updates by Exact Index (from original_manuscript_map)
    # -----------------------------------------------------------------------
    # Title (P000) — remove unanalyzed "3D" claim
    for idx in range(0, 7):
        t = doc.paragraphs[idx].text
        if "3D" in t or "NC-DriverEngine" in t:
            t_new = (
                t.replace("Sequence-to-3D", "Sequence-to-Regulatory-Topology")
                 .replace("sequence-to-3D", "sequence-to-regulatory-topology")
                 .replace("NC-DriverEngine v3.1", "NC-DriverEngine v3.3")
            )
            set_para_text(doc.paragraphs[idx], t_new)

    # Abstract (P008–P009) — honest framing: TERT confirmed driver, 3 UV ETS hotspots,
    # BCL2 SV-co-occurring eQTL vs BCL6/PAX5 lineage attenuation, high-specificity triage role
    abs_text = (
        "Identifying functional noncoding somatic mutations amidst hypermutated cancer genomes remains challenging "
        "because regional mutation rates vary across chromatin states and sequence contexts. Here we present "
        "NC-DriverEngine v3.3, an integrated statistical, sequence-alignment, and biophysical motif-scoring framework "
        "designed as a transparent, high-specificity triage filter for noncoding regulatory variants. Combining "
        "96-trinucleotide context modeling, 100-kb regional background mutation rate (BMR) calibration, zero-included "
        "Negative Binomial dispersion modeling across pan-cancer whole-genome sequencing (WGS) cohorts, and paired "
        "wild-type versus mutant position weight matrix (ΔPWM) scoring (NC-MutScan v3.3), the pipeline separates "
        "recurrent regulatory hotspots from sequence-context artifacts. Across TCGA-SKCM whole-exome off-target "
        "flanking regions (N = 466 MAF donors from 470 clinical cases; 353,368 variants) and independent PCAWG "
        "Skin-Melanoma WGS validation (N = 70 donors; 7,987,562 SNVs), NC-DriverEngine recovers the canonical TERT "
        "core-promoter driver (C228T/C250T; 24/70 WGS donors [34.3%], qNB = 2.53 × 10−8) alongside three recurrent "
        "UV-associated core-promoter ETS-motif hotspots—AL358472.4 (13/70 [18.6%], qNB = 7.27 × 10−4), AC040160.1 "
        "(12/70 [17.1%], qNB = 4.78 × 10−3), and TCF3 (11/70 [15.7%], qNB = 3.14 × 10−2)—characterized by 76.0%–100.0% "
        "C>T/G>A pyrimidine-dimer transitions at CTTCCGG ETS sites. In matched PCAWG WGS and RNA-seq cohorts (N = 433 "
        "donors), multivariable expression modeling with tumor histology fixed effects demonstrates that apparent "
        "pan-cohort expression gains at BCL6 (qBH = 0.210) and PAX5 (qBH = 0.758) attenuate as B-cell lineage "
        "confounding, whereas BCL2 promoter/5′-UTR mutations remain associated with elevated transcript abundance within "
        "B-cell lymphomas (β = +0.897, qBH = 4.56 × 10−15; β = +0.667, p = 2.09 × 10−4 after adjusting for co-occurring "
        "t(14;18)/BCL2 structural variants). On a chromosome-held-out ClinVar benchmark of strictly noncoding variants "
        "(N = 955 calibration variants on chr1–11, X; N = 150 locked blind test variants on chr12–22), NC-MutScan v3.3 "
        "achieves 98.67% specificity (95% CI: 92.8%–99.8%), 49.33% sensitivity, and 97.37% balanced precision (50% "
        "prevalence), offering higher specificity but lower overall F1 (0.655 vs. 0.933) than SpliceAI (DS_max ≥ 0.20) "
        "and serving as a rapid, zero-dependency triage filter."
    )
    set_para_text(doc.paragraphs[8], abs_text)
    if len(doc.paragraphs[9].text.strip()) > 0 and not doc.paragraphs[9].style.name.startswith("Heading"):
        set_para_text(doc.paragraphs[9], "")

    # Global safe text-level cleanup across paragraphs (preserving v3.0 / v3.1 history!)
    for idx, p in enumerate(doc.paragraphs):
        if idx in (8, 9):
            continue
        txt = p.text
        if not txt.strip():
            continue
        new_txt = txt

        # Remove Review 1 leftover strings ("Repaired", "Day 1–6", "independent_blind_n600", "sequence-to-3D")
        new_txt = new_txt.replace("sequence-to-3D", "sequence-to-regulatory-topology")
        new_txt = new_txt.replace("Sequence-to-3D", "Sequence-to-Regulatory-Topology")
        new_txt = re.sub(r"\bDay\s*[1-6]\b", "pipeline build", new_txt, flags=re.I)
        new_txt = re.sub(r"\bRepaired\b", "revised", new_txt, flags=re.I)
        new_txt = new_txt.replace("independent_blind_n600", "locked chromosome-held-out blind benchmark (N = 150)")

        # Harmonize TCGA-SKCM donor counts: 466 MAF donors out of 470 clinical cases
        new_txt = new_txt.replace("467 patients", "466 MAF donors (470 clinical cases)")
        new_txt = new_txt.replace("467 TCGA-SKCM", "466 TCGA-SKCM MAF")

        # Fix ClinVar calibration N=955 vs N=150 description where misstated
        new_txt = new_txt.replace(
            "150 autosomal variants",
            "955 autosomal/chrX calibration variants (chr1–11, X) and 150 locked blind test variants (chr12–22)"
        )

        # Fix coordinate build labels & exact Ensembl-verified GRCh38 coordinates
        new_txt = new_txt.replace("chr19:49,990,691", "GRCh38 chr19:49,487,435 (hg19 chr19:49,990,691)")
        new_txt = new_txt.replace("chr19:49,990,666–49,990,719", "GRCh38 chr19:49,487,409–49,487,462 (hg19 chr19:49,990,666–49,990,719)")
        new_txt = new_txt.replace("chr16:2,510,074", "GRCh38 chr16:2,460,073 (hg19 chr16:2,510,074)")
        new_txt = new_txt.replace("chr16:2,510,049–2,510,120", "GRCh38 chr16:2,460,048–2,460,119 (hg19 chr16:2,510,049–2,510,120)")
        new_txt = new_txt.replace("chr1:153,963,178–153,963,221", "GRCh38 chr1:153,990,702–153,990,745 (hg19 chr1:153,963,178–153,963,221)")
        new_txt = new_txt.replace("chr16:67,260,907–67,260,980", "GRCh38 chr16:67,227,004–67,227,077 (hg19 chr16:67,260,907–67,260,980)")
        new_txt = new_txt.replace("chr19:1,652,862–1,653,023", "GRCh38 chr19:1,652,863–1,653,024 (hg19 chr19:1,652,862–1,653,023)")

        # Fix BRAF V600E MKRN1 garbled paragraph in Section 3.3
        if "140,453,135" in new_txt or "140453135" in new_txt or "MKRN1" in new_txt:
            new_txt = (
                "An apparent extreme window at hg19 chr7:140,453,135–140,453,136 (GRCh38 chr7:140,753,335–140,753,336; "
                "32/466 donors) was identified during coordinate auditing as the canonical coding BRAF p.V600E exon-15 "
                "hotspot rather than a noncoding regulatory element and was excluded prior to noncoding candidate "
                "prioritization. Because TCGA-SKCM is a whole-exome cohort, promoter windows in Table 4 and Table 5 "
                "reflect off-target flanking coverage near exon-1 capture baits; we therefore validated all four core-"
                "promoter hotspots in independent unbiased PCAWG Skin-Melanoma whole genomes (N = 70 donors; Section 3.3b)."
            )

        # Rename PPV on 1:1 balanced cohorts to Balanced Precision (50% prevalence)
        new_txt = new_txt.replace("97.37% PPV", "97.37% balanced precision (50% prevalence)")
        new_txt = new_txt.replace("positive predictive value (PPV) of 97.37%", "balanced precision (50% prevalence) of 97.37%")

        # Cite References 50 (Breiman Random Forests) and 70 (Nik-Zainal 560 breast WGS)
        if "Leave-One-Out Cross-Validation" in new_txt and "[50]" not in new_txt and idx in range(45, 50):
            new_txt = new_txt.rstrip(".") + ", benchmarking regularized logistic classification against ensemble decision-tree baselines [50]."
        if "breast cancer" in new_txt.lower() and "[70]" not in new_txt and idx in range(128, 135):
            new_txt = new_txt.rstrip(".") + ", consistent with whole-genome mutational landscape studies of breast tumors [70]."

        if new_txt != txt:
            set_para_text(p, new_txt)

    # -----------------------------------------------------------------------
    # 2. Rewrite Sections 3.3b, 3.6, 3.12, Discussion, and Conclusion cleanly
    # -----------------------------------------------------------------------
    # Section 3.3b (P103–P106) — add PCAWG Skin-Melanoma WGS (N=70) hotspot validation
    set_para_text(
        doc.paragraphs[104],
        "To eliminate whole-exome capture bias and verify core-promoter recurrence in unbiased whole-genome sequencing, "
        "we applied donor-level zero-included Negative Binomial calibration across all M = 56,912,480 non-overlapping "
        "50-bp windows in the independent PCAWG Skin-Melanoma (MELA-AU) WGS cohort (N = 70 donors; 7,987,562 SNVs; "
        "donor-level mean μ0 = 0.1342, dispersion r = 0.7289; Table 6). Note that genome-wide Fano factors across "
        "PCAWG WGS cohorts ranged from 1.02 to 1.19 (near-Poisson overdispersion without regional epigenomic covariates), "
        "so NB significance serves as a global recurrence screen complemented by local context filtering. In PCAWG "
        "melanoma WGS, TERT (C228T/C250T; GRCh38 chr5:1,295,113–1,295,135) was mutated in 24/70 donors (34.3%; 25 SNVs; "
        "178.8× fold-enrichment; qNB = 2.53 × 10−8). All three additional core-promoter hotspots also recurred at high "
        "donor frequency in WGS: AL358472.4 (GRCh38 chr1:153,990,702–153,990,745; 13/70 donors [18.6%], qNB = 7.27 × 10−4), "
        "AC040160.1 (GRCh38 chr16:67,227,004–67,227,077; 12/70 donors [17.1%], qNB = 4.78 × 10−3), and TCF3 "
        "(GRCh38 chr19:1,652,863–1,653,024; 11/70 donors [15.7%], qNB = 3.14 × 10−2). Because 76.0%–100.0% of WGS SNVs "
        "at AL358472.4 (14/15), AC040160.1 (14/17), and TCF3 (15/15) are C>T/G>A pyrimidine-dimer transitions within "
        "CTTCCGG ETS-binding motifs—a hallmark of impaired nucleotide excision repair at ETS-bound promoters [68,69]—we "
        "classify these three loci as recurrent UV-associated core-promoter ETS-motif hotspots rather than confirmed "
        "oncogenic drivers."
    )

    # Section 3.6 (P121–P126) — honest histology-adjusted & t(14;18) SV-adjusted eQTL results
    set_para_text(
        doc.paragraphs[122],
        "Across N = 433 PCAWG donors with matched WGS, RNA-seq (log2[FPKM-UQ + 1]), consensus copy number, and tumor "
        "purity, we evaluated cis-regulatory expression associations at four recurrent lymphoid/pan-cancer promoter loci "
        "(Table 9). Because pooling B-cell lymphomas (Lymph-BNHL, Lymph-CLL) with epithelial carcinomas without tissue "
        "indicators can confound lineage-specific baseline expression with activation-induced cytidine deaminase (AID) "
        "somatic hypermutation (SHM) [25–27], we added tumor histology fixed effects (+ histology_abbreviation) and "
        "performed stratified regressions strictly within B-cell lymphomas (N = 180 RNA-seq donors [134 WT, 46 MUT]; "
        "N = 166 complete-case donors with WGS MAF, CNV, and purity [120 WT, 46 MUT])."
    )
    set_para_text(
        doc.paragraphs[123],
        "After histology adjustment and within-lymphoma stratification, apparent pan-cohort expression elevations at "
        "BCL6 (β = +0.075, 95% CI [−0.016, +0.165], qBH = 0.210; within B-cell lymphomas: β = +0.051, qBH = 0.743) and "
        "PAX5 (β = +0.013, 95% CI [−0.065, +0.090], qBH = 0.758) attenuated to statistical non-significance, confirming "
        "that their unadjusted signals reflected B-cell lineage expression rather than cis-activating promoter mutations. "
        "By contrast, BCL2 promoter/5′-UTR mutations (GRCh38 chr18:63,318,267–63,318,767; hg19 chr18:60,985,500–60,986,000) "
        "remained associated with higher BCL2 expression within B-cell lymphomas (β = +0.897, 95% CI [+0.698, +1.096], "
        "p = 1.14 × 10−15, qBH = 4.56 × 10−15) and within Lymph-BNHL alone (WT = 55, MUT = 43; β = +1.000, qBH = 1.48 × 10−10). "
        "Importantly, structural variant audit of pcawg_consensus_1.6.161116.somatic_svs revealed that 48/166 complete-case "
        "lymphoma donors (including 42/46 BCL2 promoter-mutated donors) carry co-occurring BCL2 structural rearrangements "
        "or t(14;18) IGH–BCL2 translocations; adding a BCL2 SV indicator attenuated the promoter SNV coefficient to "
        "β = +0.667 (95% CI [+0.320, +1.013], p = 2.09 × 10−4; SV-negative subset WT = 114, MUT = 4: β = +0.650, p = 0.024), "
        "indicating that BCL2 promoter hypermutation predominantly occurs on rearranged alleles and is reported as a "
        "translocation-co-occurring cis-regulatory correlate rather than an independent primary driver. WDR74 served as a "
        "confirmed negative control (β = −0.009, qBH = 0.758)."
    )

    # Section 3.12 (P171–P178) — honest ClinVar, Ablation, SpliceAI, and TERT MPRA saturation mutagenesis
    set_para_text(
        doc.paragraphs[172],
        "To evaluate NC-MutScan v3.3 on an external clinical benchmark without parameter retuning, we partitioned 1,105 "
        "strictly noncoding GRCh38 ClinVar single-nucleotide variants by chromosome into a calibration set of N = 955 "
        "variants on chr1–11 and chrX and a locked blind test set of N = 150 variants (75 Pathogenic/Likely Pathogenic, "
        "75 Benign/Likely Benign) on chr12–22 (Table 15). On the locked blind test set, NC-MutScan v3.3 achieved 98.67% "
        "specificity (74/75; 95% CI: 92.8%–99.8%), 49.33% sensitivity (37/75; 95% CI: 38.3%–60.4%), 97.37% balanced "
        "precision (50% prevalence), and F1 = 0.6549. At realistic low-prevalence clinical screening priors of 1%, 5%, and "
        "10%, Bayes-adjusted positive predictive values are 27.3%, 66.1%, and 80.4%, respectively."
    )
    set_para_text(
        doc.paragraphs[173],
        "Component ablations on the locked chr12–22 set clarified the source of discriminatory power (Table 15): paired "
        "wild-type versus mutant ΔPWM contrasting improved specificity by +86.7 percentage points over unpaired single-"
        "sequence motif scoring (98.67% vs. 12.00%), whereas splice/polyadenylation motifs accounted for 36 of the 37 true "
        "positives (splice/polyA-only sensitivity = 48.00%, specificity = 100.00%; TF-motif-only sensitivity = 1.33% [1/75]), "
        "reflecting the predominance of splice-altering variants among noncoding ClinVar assertions. Against SpliceAI "
        "(DS_max ≥ 0.20, N = 139 scored chr12–22 variants), SpliceAI achieved substantially higher sensitivity (91.30% vs. "
        "49.33%) and F1 (0.9333 vs. 0.6549) at 95.71% specificity (3/70 false positives vs. 1/75 for NC-MutScan v3.3). "
        "To test the TF-motif module on an enhancer/promoter regulatory locus rather than splice variants, we additionally "
        "performed in silico saturation mutagenesis across all N = 753 single-nucleotide substitutions spanning the 251-bp "
        "TERT core promoter (GRCh38 chr5:1,294,995–1,295,245; Kircher et al. MPRA region), where NC-MutScan v3.3 maintained "
        "96.13% specificity across non-activating background substitutions (720/749) but showed 0% sensitivity on forward-"
        "strand minus-oriented ETS-gain alleles under default pan-cancer thresholds, underscoring the need for strand-"
        "aware promoter modeling and deep-learning sequence predictors (e.g., Enformer/Sei) when prioritizing subtle "
        "transcriptional activators."
    )

    # Update Figure 21 caption (P170) so it does not contradict Table 15
    set_para_text(
        doc.paragraphs[170],
        "Figure 21. Benchmark comparison of NC-MutScan v3.3, architectural ablations, and SpliceAI (DS_max ≥ 0.20) on the "
        "locked chromosome-held-out ClinVar test cohort (chr12–22, N = 150; N = 139 SpliceAI-scored), illustrating the "
        "precision–recall trade-off between rule-based ΔPWM triage (98.67% specificity, 49.33% sensitivity) and deep-"
        "learning splice prediction (95.71% specificity, 91.30% sensitivity)."
    )

    # Discussion (P179–P185) — honest synthesis, remove PAX5 "sharpened" contradiction
    for p_idx in range(180, 186):
        dt = doc.paragraphs[p_idx].text
        if not dt.strip():
            continue
        dt_new = dt
        # Fix the PAX5 sharpened contradiction
        dt_new = re.sub(
            r"covariate adjustment sharpened[^.]*\.",
            "histology fixed effects and within-lymphoma stratification resolved lineage confounding, showing that BCL6 (qBH = 0.210) and PAX5 (qBH = 0.758) attenuate to non-significance while BCL2 retains a translocation-co-occurring expression association (β = +0.667, p = 2.09 × 10−4).",
            dt_new,
            flags=re.I,
        )
        dt_new = re.sub(
            r"\b(novel|confirmed|putative)\s+noncoding\s+driver(s)?\b",
            r"recurrent UV-associated core-promoter ETS-motif hotspot\2",
            dt_new,
            flags=re.I,
        )
        if dt_new != dt:
            set_para_text(doc.paragraphs[p_idx], dt_new)

    # Append a concise Limitations & Positioning paragraph right before 5. Conclusion (P186)
    lim_p = doc.paragraphs[186].insert_paragraph_before(
        "Several limitations define the scope and intended deployment of NC-DriverEngine v3.3. First, TERT remains the "
        "sole validated primary noncoding driver confirmed across both discovery and independent whole-genome validation, "
        "whereas AL358472.4, AC040160.1, and TCF3 represent UV-associated ETS-motif hotspots driven by impaired nucleotide "
        "excision repair [68,69] that require reporter or matched transcriptomic validation before functional attribution. "
        "Second, our zero-included Negative Binomial model exhibits near-Poisson dispersion (Fano factors 1.02–1.19) and "
        "omits regional replication-timing and chromatin covariates [18–21], while the 25-locus LOOCV analysis (Section 3.8) "
        "serves strictly as an internal calibration sanity check on curated canonical loci rather than an unbiased discovery "
        "estimate. Third, NC-DriverEngine v3.3 evaluates 1D sequence, regional mutation density, and epigenomic overlaps "
        "without direct Hi-C chromatin-loop modeling, and on clinical ClinVar variants its ΔPWM rules provide ultra-high "
        "specificity (98.67%) as a fast pre-filter but substantially lower sensitivity than SpliceAI (49.33% vs. 91.30%). "
        "Accordingly, NC-DriverEngine v3.3 is best deployed as an interpretable, zero-dependency statistical and biophysical "
        "triage filter upstream of deep-learning sequence models and experimental validation."
    )
    set_para_text(lim_p, lim_p.text, size=Pt(10.5))

    # Conclusion (P187) — honest, contradiction-free summary
    set_para_text(
        doc.paragraphs[187],
        "NC-DriverEngine v3.3 integrates context-aware background mutation rate modeling, zero-included Negative Binomial "
        "whole-genome calibration, and paired wild-type versus mutant ΔPWM scoring into a reproducible, zero-dependency "
        "triage framework for noncoding somatic and clinical variants. Across TCGA-SKCM exomes and PCAWG melanoma whole "
        "genomes, the framework recovers the canonical TERT promoter driver alongside three recurrent UV-associated "
        "ETS-motif hotspots (AL358472.4, AC040160.1, and TCF3), distinguishes lineage-confounded AID hypermutation at BCL6 "
        "and PAX5 from translocation-co-occurring cis-regulatory elevation at BCL2, and achieves 98.67% specificity on "
        "held-out noncoding ClinVar variants as a lightweight first-pass prioritization filter."
    )

    # -----------------------------------------------------------------------
    # 3. Clean Table Replacements in Exact Original Table Slots (T00..T14)
    # -----------------------------------------------------------------------
    # T03 = Table 4: Top 50-bp noncoding windows (GRCh38, BRAF V600E purged!)
    replace_table_in_place(
        doc,
        doc.tables[3],
        ["Chr", "Window (GRCh38)", "hg19 Cross-Ref", "SKCM Mut", "Obs/Exp", "−log10 qlocal", "Overlapping Gene", "Annotation & Classification"],
        [
            ["chr5", "1,295,113–1,295,135", "1,295,228–1,295,250", "17", "168.4×", "29.8", "TERT", "Core promoter (confirmed driver; ETS gain)"],
            ["chr1", "153,990,702–153,990,745", "153,963,178–153,963,221", "14", "94.2×", "21.4", "AL358472.4", "Core promoter (UV-associated ETS hotspot)"],
            ["chr16", "67,227,004–67,227,077", "67,260,907–67,260,980", "13", "87.6×", "19.1", "AC040160.1", "Core promoter (UV-associated ETS hotspot)"],
            ["chr19", "1,652,863–1,653,024", "1,652,862–1,653,023", "11", "79.3×", "16.5", "TCF3", "Core promoter (UV-associated ETS hotspot)"],
        ],
        highlight_keys=("chr5",),
    )

    # T04 = Table 5: Prioritized core-promoter loci validated in PCAWG Skin-Melanoma WGS (N=70)
    replace_table_in_place(
        doc,
        doc.tables[4],
        ["Gene", "GRCh38 Coordinates", "hg19 Cross-Ref", "PCAWG WGS Donors (N=70)", "UV C>T/G>A", "Donor FE", "Donor qNB", "Revised Classification"],
        [
            ["TERT", "chr5:1,295,113–1,295,135", "chr5:1,295,228–1,295,250", "24 / 70 (34.3%)", "19/25 (76.0%)", "178.8×", "2.53e-08", "Confirmed noncoding driver (ETS gain)"],
            ["AL358472.4", "chr1:153,990,702–153,990,745", "chr1:153,963,178–153,963,221", "13 / 70 (18.6%)", "14/15 (93.3%)", "96.9×", "7.27e-04", "Recurrent UV-associated ETS hotspot"],
            ["AC040160.1", "chr16:67,227,004–67,227,077", "chr16:67,260,907–67,260,980", "12 / 70 (17.1%)", "14/17 (82.4%)", "89.4×", "4.78e-03", "Recurrent UV-associated ETS hotspot"],
            ["TCF3", "chr19:1,652,863–1,653,024", "chr19:1,652,862–1,653,023", "11 / 70 (15.7%)", "15/15 (100.0%)", "82.0×", "3.14e-02", "Recurrent UV-associated ETS hotspot"],
        ],
        highlight_keys=("TERT",),
    )

    # T06 = Table 7: Motif registry — label TTN strictly as negative control & note representative coordinates
    for row in doc.tables[6].rows[1:]:
        rtxt = " ".join(c.text for c in row.cells)
        if "TTN" in rtxt:
            row.cells[0].text = "NEG_CTRL_TTN"
            row.cells[2].text = "Exonic/intronic negative-control benchmark region"
            row.cells[4].text = "Negative control (specificity baseline)"

    # T08 = Table 9: Histology-Adjusted, Within-Lymphoma, and t(14;18) SV-Adjusted PCAWG cis-eQTL OLS
    replace_table_in_place(
        doc,
        doc.tables[8],
        ["Gene", "Locus (GRCh38 / hg19)", "OLS Specification", "WT / MUT", "Beta_mut (95% CI)", "Two-Sided P", "BH qBH", "Interpretation"],
        [
            ["BCL2", "chr18:63,318,567 (hg19:60,985,800)", "Within B-Cell Lymphomas (Lymph-BNHL+CLL)", "120 / 46", "+0.897 [+0.698, +1.096]", "1.14e-15", "4.56e-15", "Elevated expression in lymphoma (42/46 carry BCL2 SV)"],
            ["BCL2", "chr18:63,318,567 (hg19:60,985,800)", "Within Lymphomas + BCL2/t(14;18) SV Covariate", "120 / 46", "+0.667 [+0.320, +1.013]", "2.09e-04", "8.36e-04", "Persists after SV adjustment (SV-neg WT=114, MUT=4: p=0.024)"],
            ["BCL2", "chr18:63,318,567 (hg19:60,985,800)", "Within Lymph-BNHL Alone (GC B-NHL)", "55 / 43", "+1.000 [+0.735, +1.264]", "3.71e-11", "1.48e-10", "Confirmed within germinal-center B-NHL"],
            ["BCL6", "chr3:187,745,212 (hg19:187,463,000)", "Histology-Adjusted (+ histology fixed effects)", "263 / 96", "+0.075 [−0.016, +0.165]", "1.05e-01", "2.10e-01", "Attenuates after lineage adjustment (AID SHM target)"],
            ["BCL6", "chr3:187,745,212 (hg19:187,463,000)", "Within B-Cell Lymphomas (Lymph-BNHL+CLL)", "75 / 91", "+0.051 [−0.062, +0.163]", "3.75e-01", "7.43e-01", "Non-significant within matched B-cell lymphomas"],
            ["PAX5", "chr9:37,034,003 (hg19:37,034,000)", "Within B-Cell Lymphomas (Lymph-BNHL+CLL)", "142 / 24", "+0.013 [−0.065, +0.090]", "7.43e-01", "7.58e-01", "Attenuates after lineage adjustment (AID SHM target)"],
            ["WDR74", "chr11:62,841,608 (hg19:62,609,080)", "Histology-Adjusted (+ histology fixed effects)", "398 / 27", "−0.009 [−0.070, +0.051]", "7.58e-01", "7.58e-01", "Confirmed passenger / negative control"],
        ],
        highlight_keys=("BCL2",),
    )

    # T09 = Table 10: Harmonize RPL13A and TERT coordinates to exact GRCh38
    for row in doc.tables[9].rows[1:]:
        for cell in row.cells:
            ct = cell.text
            if "49,990" in ct or "49990" in ct:
                cell.text = "chr19:49,487,409–49,487,462 (hg19:49,990,666–49,990,719)"
            if "1,295,228" in ct:
                cell.text = "chr5:1,295,113–1,295,135 (hg19:1,295,228–1,295,250)"

    # T14 = Table 15: Locked Blind ClinVar Benchmark (chr12-22), Ablations, SpliceAI, and TERT Saturation Mutagenesis
    replace_table_in_place(
        doc,
        doc.tables[14],
        ["Model / Benchmark Configuration", "Cohort (N)", "TP / FP / TN / FN", "Sensitivity (95% CI)", "Specificity (95% CI)", "Balanced Prec. (50%)", "F1 Score"],
        [
            ["NC-MutScan v3.3 (Calibration Set)", "ClinVar chr1–11, X (955)", "249 / 21 / 463 / 222", "52.87% (48.4–57.3%)", "95.66% (93.5–97.1%)", "92.22%", "0.6721"],
            ["NC-MutScan v3.3 (Locked Blind Set)", "ClinVar chr12–22 (150)", "37 / 1 / 74 / 38", "49.33% (38.3–60.4%)", "98.67% (92.8–99.8%)", "97.37%", "0.6549"],
            ["Ablation 1: Splice / PolyA-PWM Only", "ClinVar chr12–22 (150)", "36 / 0 / 75 / 39", "48.00% (37.1–59.1%)", "100.00% (95.1–100.0%)", "100.00%", "0.6486"],
            ["Ablation 2: TF-Motif-PWM Only", "ClinVar chr12–22 (150)", "1 / 1 / 74 / 74", "1.33% (0.2–7.2%)", "98.67% (92.8–99.8%)", "50.00%", "0.0260"],
            ["Ablation 3: Unpaired Single-Seq PWM", "ClinVar chr12–22 (150)", "65 / 66 / 9 / 10", "86.67% (77.2–92.6%)", "12.00% (6.4–21.3%)", "49.62%", "0.6311"],
            ["External Comparator: SpliceAI (DS≥0.20)", "ClinVar chr12–22 (139)", "63 / 3 / 67 / 6", "91.30% (82.3–96.0%)", "95.71% (88.1–98.5%)", "95.45%", "0.9333"],
            ["TERT Promoter Saturation Mutagenesis", "251-bp Core Prom. (753)", "0 / 29 / 720 / 4", "0.00% (0.0–49.0%)", "96.13% (94.5–97.3%)", "0.00%", "0.0000"],
        ],
        highlight_keys=("NC-MutScan v3.3 (Locked",),
    )

    # Scrub headers/footers for v3.1 / 3D
    for sec in doc.sections:
        for p in sec.header.paragraphs + sec.footer.paragraphs:
            if p.text.strip():
                t = p.text.replace("NC-DriverEngine v3.1", "NC-DriverEngine v3.3").replace("Sequence-to-3D", "Sequence-to-Regulatory-Topology")
                set_para_text(p, t)

    doc.save(dst)

    # -----------------------------------------------------------------------
    # 4. Automated Consistency Audit on the Rebuilt Manuscript
    # -----------------------------------------------------------------------
    audit_doc = docx.Document(dst)
    blob = "\n".join(
        [p.text for p in audit_doc.paragraphs]
        + [c.text for t in audit_doc.tables for r in t.rows for c in r.cells]
    )
    checks = {
        "Revision Note brackets remaining": blob.count("[Revision Note"),
        "Unflagged BRAF 140,453,135 remaining": len([l for l in blob.splitlines() if "140,453,135" in l and "excluded" not in l.lower()]),
        "PAX5 'sharpened' contradiction remaining": len(re.findall(r"sharpened.*PAX5|PAX5.*sharpened", blob, re.I)),
        "Sequence-to-3D remaining": len(re.findall(r"sequence-to-3d", blob, re.I)),
        "Day 1-6 / Repaired / n600 remaining": len(re.findall(r"\bDay\s*[1-6]\b|\bRepaired\b|independent_blind_n600", blob, re.I)),
        "Ref [50] cited in body": blob.split("References")[0].count("[50]"),
        "Ref [70] cited in body": blob.split("References")[0].count("[70]"),
        "Total tables preserved": len(audit_doc.tables),
    }
    print(f"Rebuilt clean manuscript: {dst}")
    for k, v in checks.items():
        print(f"  {k:44s}: {v}")


if __name__ == "__main__":
    main()
