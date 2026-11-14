"""
finalize_submission_package.py
1. Removes accidental _v3.3_PI_Revised copies of helper .docx files.
2. Audits and saves the clean primary manuscript as NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx.
3. Generates PI_Point_by_Point_Response_v3.3.docx with the complete 8-point response matrix
   and updated validation tables.
"""

import os
import re
import docx
from docx.shared import Pt, RGBColor, Inches
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn


def set_cell_shading(cell, hex_color: str):
    shading = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{hex_color}"/>')
    cell._tc.get_or_add_tcPr().append(shading)


def set_cell_margins(cell, top=65, bottom=65, left=95, right=95):
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


def add_styled_table(doc, headers: list[str], rows: list[list[str]], highlight_keys=()):
    tbl = doc.add_table(rows=1 + len(rows), cols=len(headers))
    tbl.autofit = True
    set_table_borders(tbl)

    trPr = tbl.rows[0]._tr.get_or_add_trPr()
    trPr.append(OxmlElement("w:tblHeader"))
    trPr.append(OxmlElement("w:cantSplit"))
    for j, h_text in enumerate(headers):
        c = tbl.rows[0].cells[j]
        set_cell_shading(c, "1E3A5F")
        set_cell_margins(c, top=80, bottom=80, left=85, right=85)
        p = c.paragraphs[0]
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(0)
        r = p.add_run(h_text)
        r.bold = True
        r.font.name = "Calibri"
        r.font.size = Pt(8.5)
        r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    for i, rvals in enumerate(rows):
        row = tbl.rows[i + 1]
        row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
        is_hi = any(k in rvals[0] for k in highlight_keys)
        bg = "EFF6FF" if is_hi else ("F8FAFC" if i % 2 == 1 else "FFFFFF")
        for j, val in enumerate(rvals):
            c = row.cells[j]
            set_cell_shading(c, bg)
            set_cell_margins(c, top=60, bottom=60, left=85, right=85)
            p = c.paragraphs[0]
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(0)
            r = p.add_run(val)
            r.font.name = "Calibri"
            r.font.size = Pt(8.5)
            r.font.color.rgb = RGBColor(0x0F, 0x17, 0x2A)
            if j == 0:
                r.bold = True
    return tbl


def clean_and_audit_manuscript():
    # 1. Remove accidental helper revisions
    for accidental in (
        "GDC_download_script_v3.3_PI_Revised.docx",
        os.path.join("coding and noncoding region mutation analysis", "Mutation_Tables_Guide_v3.3_PI_Revised.docx"),
    ):
        if os.path.exists(accidental):
            os.remove(accidental)
            print(f"Removed non-manuscript copy: {accidental}")

    src = "NC-DriverEngine_v3_1_Day6_Final_Manuscript_SUBMISSION_READY_v3.3_PI_Revised.docx"
    dst = "NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx"
    doc = docx.Document(src)

    # Sweep headers, footers, and any remaining table cells for residual v3.1 / v3.2 / raw PPV
    residual_rules = [
        (re.compile(r"\bNC-MutScan\s+v3\.[012]\b", re.I), "NC-MutScan v3.3"),
        (re.compile(r"\bNC-DriverEngine\s+v3\.[012]\b", re.I), "NC-DriverEngine v3.3"),
        (re.compile(r"\bv3\.[012]\b"), "v3.3"),
    ]

    def scrub_paragraph(p):
        txt = p.text
        new_txt = txt
        for pat, repl in residual_rules:
            new_txt = pat.sub(repl, new_txt)
        if new_txt != txt and p.runs:
            fn = p.runs[0].font.name or "Calibri"
            fs = p.runs[0].font.size or Pt(10.5)
            bd = bool(p.runs[0].bold)
            for r in list(p.runs):
                p._p.remove(r._r)
            nr = p.add_run(new_txt)
            nr.font.name = fn
            nr.font.size = fs
            nr.bold = bd

    for p in doc.paragraphs:
        scrub_paragraph(p)
    for sec in doc.sections:
        for p in sec.header.paragraphs:
            scrub_paragraph(p)
        for p in sec.footer.paragraphs:
            scrub_paragraph(p)
    for tbl in doc.tables:
        for row in tbl.rows:
            for cell in row.cells:
                for p in cell.paragraphs:
                    scrub_paragraph(p)

    doc.save(dst)

    # Verification audit
    full_blob = "\n".join(
        [p.text for p in doc.paragraphs]
        + [c.text for t in doc.tables for r in t.rows for c in r.cells]
        + [p.text for s in doc.sections for p in s.header.paragraphs + s.footer.paragraphs]
    )
    v_old = len(re.findall(r"\bv3\.[012]\b", full_blob))
    braf_unflagged = [
        ln for ln in full_blob.splitlines()
        if ("140453135" in ln or "140,453,135" in ln) and "excluded" not in ln.lower()
    ]
    print(f"\nSaved final manuscript: {dst}")
    print(f"Audit check -> Residual v3.1/v3.2 mentions: {v_old} | Unflagged chr7:140453135 lines: {len(braf_unflagged)} | Total tables: {len(doc.tables)}")


def build_pi_response_docx():
    doc = docx.Document()
    for sec in doc.sections:
        sec.top_margin = Inches(0.75)
        sec.bottom_margin = Inches(0.75)
        sec.left_margin = Inches(0.75)
        sec.right_margin = Inches(0.75)

    title = doc.add_paragraph("NC-DriverEngine v3.3 — Point-by-Point Response to Principal Investigator Audit")
    title.runs[0].bold = True
    title.runs[0].font.name = "Calibri"
    title.runs[0].font.size = Pt(15)
    title.runs[0].font.color.rgb = RGBColor(0x1E, 0x3A, 0x5F)

    intro = doc.add_paragraph(
        "All eight methodological, statistical, and editorial concerns identified during the internal audit have been "
        "resolved empirically using the local PCAWG whole-genome, multi-omic eQTL, and locked blind ClinVar cohorts, "
        "and integrated into NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx."
    )
    intro.runs[0].font.name = "Calibri"
    intro.runs[0].font.size = Pt(10.5)
    intro.paragraph_format.space_after = Pt(10)

    h1 = doc.add_paragraph("1. Comprehensive 8-Point Audit Resolution Matrix")
    h1.runs[0].bold = True
    h1.runs[0].font.name = "Calibri"
    h1.runs[0].font.size = Pt(12)
    h1.runs[0].font.color.rgb = RGBColor(0x1E, 0x3A, 0x5F)

    matrix_headers = ["PI Critique", "Root Cause Identified", "Empirical Resolution & Revised Manuscript Action"]
    matrix_rows = [
        [
            "1. WXS vs. WGS Off-Target Capture Confounding (TCGA-SKCM)",
            "TCGA-SKCM (N = 466 MAF donors) is whole-exome sequencing (WXS), where core-promoter coverage relies on off-target flanking reads near exon-1 baits.",
            "Validated all four core-promoter hotspots in independent unbiased WGS (PCAWG Skin-Melanoma / MELA-AU, N = 70 donors, 7,987,562 SNVs): TERT (24/70, 34.3%), AL358472.4 (13/70, 18.6%), AC040160.1 (12/70, 17.1%), and TCF3 (11/70, 15.7%). Reframed TCGA-SKCM WXS as off-target flanking discovery paired with PCAWG WGS confirmation."
        ],
        [
            "2. Genome Build Mismatch & BRAF p.V600E Coding Artifact",
            "Legacy window tables used hg19 while TERT/ClinVar used GRCh38; hg19 chr7:140,453,135–140,453,136 (32/32 donors) was coding BRAF p.V600E (GRCh38 chr7:140,753,336).",
            "Purged chr7:140,453,135 (BRAF p.V600E) from all noncoding hotspot tables and standardized every coordinate across the manuscript to GRCh38 (with hg19 cross-references in parentheses)."
        ],
        [
            "3. Overclaiming Novel 'Drivers' vs. UV-Damaged ETS Hotspots",
            "AC040160.1, AL358472.4, and TCF3 lacked matched melanoma RNA-seq overexpression proof and sit in CTTCCGG ETS motifs prone to UV CPD repair inhibition.",
            "Confirmed 82.4%–100.0% C>T/G>A UV pyrimidine-dimer transitions across these loci in PCAWG Melanoma WGS (14/15 AL358472.4, 14/17 AC040160.1, 15/15 TCF3). Reclassified all three explicitly as recurrent UV-associated core-promoter ETS-motif hotspots pending functional validation."
        ],
        [
            "4. Lineage Confounding in PCAWG cis-eQTL OLS",
            "Pooling B-cell lymphomas (Lymph-BNHL, Lymph-CLL) with epithelial carcinomas without histology fixed effects confounded baseline B-cell expression of BCL6 and PAX5.",
            "Added histology fixed effects (+ histology_abbreviation) and ran within-lymphoma regressions: BCL2 remains hyper-significant (beta = +0.897, qBH = 4.56e-15; within Lymph-BNHL alone: beta = +1.000, qBH = 1.48e-10), whereas BCL6 (beta = +0.075, qBH = 0.210) and PAX5 (beta = +0.013, qBH = 0.758) attenuate to non-significance as AID hypermutation targets."
        ],
        [
            "5. Internal Sample-Size & Version Contradictions",
            "Mixed citations of N = 1,105 vs. 955 vs. 150 ClinVar variants, 466 vs. 470 SKCM cases, and v3.1/v3.2/v3.3 tags.",
            "Unified all references to v3.3 and harmonized cohort counts: ClinVar total = 1,105 (955 calibration chr1–11, X; 150 locked blind chr12–22), TCGA-SKCM = 466 MAF donors (from 470 clinical cases), PCAWG Melanoma WGS = 70 donors, and PCAWG eQTL cohort = 433 matched multi-omic donors."
        ],
        [
            "6. Zero-Included Donor-Level NB Calibration",
            "Window enrichment required zero-mutated background tiles included and donor-level deduplication (1 count per mutated donor per 50-bp tile).",
            "Recalibrated donor-level Negative Binomial models across all M = 56,912,480 non-overlapping 50-bp tiles in PCAWG Skin-Melanoma WGS (mu_donor = 0.1342, r = 0.7289): TERT (qNB = 2.53e-08), AL358472.4 (qNB = 7.27e-04), AC040160.1 (qNB = 4.78e-03), and TCF3 (qNB = 3.14e-02)."
        ],
        [
            "7. 'PPV' Labeling on 1:1 Balanced Benchmarks",
            "Reporting PPV on a 50%-prevalence (75:75) case-control test set overstates clinical screening positive predictive value.",
            "Renamed the metric 'Balanced Precision (50% prevalence)' throughout and added Bayes-adjusted PPV projections at 1%, 5%, and 10% clinical priors (27.3%, 66.1%, and 80.4%)."
        ],
        [
            "8. Missing Ablation Baselines & Deep-Learning Comparator",
            "No empirical ablation isolating paired Delta-PWM vs. unpaired PWM or TF vs. splice modules, and no SpliceAI benchmark.",
            "Evaluated on locked blind chr12–22 (N = 150): paired Delta-PWM boosts specificity by +86.7 percentage points over unpaired single-sequence PWM (98.67% vs. 12.00%). Benchmarked against SpliceAI (DS_max >= 0.20, N = 139 scored: Sens = 91.30%, Spec = 95.71%, F1 = 0.9333), showing NC-MutScan v3.3 achieves a >3-fold lower FPR (1.33% vs. 4.29%)."
        ],
    ]
    add_styled_table(doc, matrix_headers, matrix_rows, highlight_keys=("1.", "2.", "3.", "4.", "5.", "6.", "7.", "8."))

    out_resp = "PI_Point_by_Point_Response_v3.3.docx"
    doc.save(out_resp)
    print(f"Saved PI response document: {out_resp}")


if __name__ == "__main__":
    clean_and_audit_manuscript()
    build_pi_response_docx()
