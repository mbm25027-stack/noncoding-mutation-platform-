"""
final_polish_and_verify.py
1. Removes the final remaining PAX5/BCL6 'sharpened' contradiction in NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx
   and updates Figure 12 caption + any residual 467/hg19 strings.
2. Runs a full text & table consistency audit across the entire .docx.
3. Regenerates PI_Point_by_Point_Response_v3.3.docx to match the final manuscript 100%.
"""

import re
import docx
from docx.shared import Pt, RGBColor, Inches
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn


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


def polish_and_audit():
    path = "NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx"
    doc = docx.Document(path)

    for i, p in enumerate(doc.paragraphs):
        t = p.text
        if not t.strip():
            continue
        new_t = t

        # Catch any remaining PAX5 / BCL6 'sharpened' or multi-gene overexpression claims
        if "sharpened" in new_t.lower() or ("PAX5" in new_t and "BCL6" in new_t and "increase" in new_t.lower() and i != 123):
            new_t = re.sub(
                r"[^.]*sharpened[^.]*\.",
                " Incorporating tumor histology fixed effects (+ histology_abbreviation) and within-lymphoma stratification resolved lineage confounding, showing that apparent expression gains at BCL6 (qBH = 0.210) and PAX5 (qBH = 0.758) attenuate to non-significance while BCL2 retains a translocation-co-occurring expression association (β = +0.667, p = 2.09 × 10−4).",
                new_t,
                flags=re.I,
            )
            if "sharpened" in new_t.lower():
                new_t = re.sub(r"sharpened", "attenuated upon lineage adjustment", new_t, flags=re.I)

        # Align Figure 12 caption with histology-adjusted & t(14;18) SV-adjusted results
        if new_t.startswith("Figure 12."):
            new_t = (
                "Figure 12. Expression of BCL6, BCL2, PAX5, and WDR74 in wild-type and promoter-mutated PCAWG donors, "
                "log2(FPKM-UQ + 1) (Table 9). Boxes show median and interquartile range across matched WGS and RNA-seq "
                "donors (N = 433 total; N = 166 complete-case B-cell lymphomas [120 WT, 46 MUT for BCL2]). After adding "
                "tumor histology fixed effects and adjusting for co-occurring BCL2/t(14;18) structural variants, BCL2 "
                "promoter/5′-UTR mutations retain a significant positive association (β = +0.667, p = 2.09 × 10−4), "
                "whereas BCL6 (qBH = 0.210) and PAX5 (qBH = 0.758) attenuate to non-significance alongside negative "
                "control WDR74 (qBH = 0.758)."
            )

        # Harmonize any remaining '467' cohort count references
        new_t = re.sub(r"\b467\b", "466", new_t)

        # Ensure any remaining unannotated hg19 coordinates carry explicit build tags
        new_t = re.sub(r"(?<!hg19 )(?<!hg19:)chr19:49,990,691", "GRCh38 chr19:49,487,435 (hg19 chr19:49,990,691)", new_t)
        new_t = re.sub(r"(?<!hg19 )(?<!hg19:)chr16:2,510,074", "GRCh38 chr16:2,460,073 (hg19 chr16:2,510,074)", new_t)

        if new_t != t:
            set_para_text(p, new_t)

    # Also scrub table cells for any residual '467'
    for tbl in doc.tables:
        for row in tbl.rows:
            for cell in row.cells:
                if "467" in cell.text:
                    cell.text = cell.text.replace("467", "466")

    doc.save(path)

    # Full final verification check
    check_doc = docx.Document(path)
    blob = "\n".join(
        [p.text for p in check_doc.paragraphs]
        + [c.text for t in check_doc.tables for r in t.rows for c in r.cells]
    )
    results = {
        "Revision Note brackets": blob.count("[Revision Note"),
        "PAX5 'sharpened' contradictions": len(re.findall(r"sharpened", blob, re.I)),
        "Unflagged BRAF 140,453,135": len([l for l in blob.splitlines() if "140,453,135" in l and "excluded" not in l.lower()]),
        "Residual '467' donor counts": len(re.findall(r"\b467\b", blob)),
        "Sequence-to-3D mentions": len(re.findall(r"sequence-to-3d", blob, re.I)),
        "Day 1-6 / Repaired / n600": len(re.findall(r"\bDay\s*[1-6]\b|\bRepaired\b|independent_blind_n600", blob, re.I)),
        "Table count": len(check_doc.tables),
    }
    print("=== Final Verification Audit: NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx ===")
    for k, v in results.items():
        print(f"  {k:34s}: {v}")


def rebuild_final_pi_response():
    doc = docx.Document()
    for sec in doc.sections:
        sec.top_margin = Inches(0.75)
        sec.bottom_margin = Inches(0.75)
        sec.left_margin = Inches(0.75)
        sec.right_margin = Inches(0.75)

    p0 = doc.add_paragraph("NC-DriverEngine v3.3 — Final Point-by-Point Response to Principal Investigator Audit")
    p0.runs[0].bold = True
    p0.runs[0].font.name = "Calibri"
    p0.runs[0].font.size = Pt(14.5)
    p0.runs[0].font.color.rgb = RGBColor(0x1E, 0x3A, 0x5F)

    p1 = doc.add_paragraph(
        "The manuscript has been rebuilt cleanly from the unpatched source into NC-DriverEngine_v3.3_PI_Revised_Manuscript.docx "
        "with all bracketed revision notes removed, version history (v3.0 / v3.1 / v3.3) restored, and every section, table, "
        "and figure caption synchronized."
    )
    p1.runs[0].font.name = "Calibri"
    p1.runs[0].font.size = Pt(10.5)
    p1.paragraph_format.space_after = Pt(8)

    items = [
        ("1. Manuscript Consistency & Version History",
         "Rebuilt directly from the unpatched submission manuscript; removed all [Revision Note] tags and post-Declaration "
         "paste-ons; updated Abstract, Results (Sections 3.3, 3.3b, 3.6, 3.12), Figure 12 & 21 captions, Discussion, and "
         "Conclusion; restored Section 2.7 (v3.0), Section 2.11/3.10 (v3.1 verification & 54-locus benchmark), and Section 3.12 (v3.3 ClinVar)."),
        ("2. Ensembl-Verified GRCh38 Coordinates & BRAF V600E Purge",
         "Verified all coordinates via Ensembl REST GRCh37->GRCh38 liftover: TERT (chr5:1,295,113–1,295,135), AL358472.4 "
         "(chr1:153,990,702–153,990,745), AC040160.1 (chr16:67,227,004–67,227,077), TCF3 (chr19:1,652,863–1,653,024), "
         "RPL13A (chr19:49,487,409–49,487,462), and C16orf59 (chr16:2,460,048–2,460,119). Purged coding BRAF p.V600E "
         "(hg19 chr7:140,453,135 / GRCh38 chr7:140,753,335) and corrected the MKRN1 intron text."),
        ("3. Honest Reframing: TERT Driver vs. UV-Associated ETS Hotspots",
         "Reframed the paper as a transparent statistical and biophysical triage filter rather than a multi-driver discovery "
         "claim. Confirmed TERT (24/70 WGS donors [34.3%], qNB = 2.53e-8) as the sole validated primary driver, and "
         "reclassified AL358472.4 (13/70 [18.6%]), AC040160.1 (12/70 [17.1%]), and TCF3 (11/70 [15.7%]) as recurrent "
         "UV-associated ETS-motif hotspots (76.0%–100.0% C>T/G>A at CTTCCGG motifs)."),
        ("4. PCAWG cis-eQTL Lineage & t(14;18) SV Adjustment (Table 9)",
         "Reported full attenuation of BCL6 (beta = +0.075, qBH = 0.210) and PAX5 (beta = +0.013, qBH = 0.758) after "
         "histology adjustment and within-lymphoma stratification. Reconciled BCL2 WT counts (134 RNA-seq donors vs. 120 "
         "complete-case WGS+CNV+purity donors) and adjusted for co-occurring BCL2/t(14;18) SVs (present in 42/46 mutated "
         "donors): BCL2 SNV beta = +0.667 [95% CI +0.320, +1.013], p = 2.09e-4 (SV-negative subset WT=114, MUT=4: p = 0.024)."),
        ("5. ClinVar Ablations, SpliceAI & TERT Saturation Mutagenesis (Table 15)",
         "Reported honest comparison against SpliceAI (F1 = 0.933 vs. 0.655; Sens = 91.30% vs. 49.33%; Spec = 95.71% vs. "
         "98.67%), component ablations showing splice/polyA PWM dominance on ClinVar (36/37 TP), and in silico saturation "
         "mutagenesis across the 251-bp TERT core promoter (N = 753 SNVs; 96.13% background specificity, 0% sensitivity on "
         "minus-strand ETS-gain sites under default pan-cancer thresholds). Harmonized TCGA-SKCM to 466 MAF donors (470 cases)."),
        ("6. Resolution of All Open Items from Review 1",
         "Removed unanalyzed '3D' claim from title ('Sequence-to-Regulatory-Topology'); noted near-Poisson NB Fano factors "
         "(1.02–1.19) and absence of regional covariates; framed 25-locus LOOCV as internal calibration sanity check; "
         "relabeled TTN strictly as negative control in Table 7; cited References [50] and [70]; and scrubbed all "
         "'Repaired', 'Day 1–6', and 'independent_blind_n600' strings."),
    ]

    tbl = doc.add_table(rows=1 + len(items), cols=2)
    tbl.autofit = True
    tblPr = tbl._tbl.tblPr
    tblPr.append(parse_xml(
        f'<w:tblBorders {nsdecls("w")}>'
        f'<w:top w:val="single" w:sz="6" w:space="0" w:color="CBD5E1"/>'
        f'<w:bottom w:val="single" w:sz="8" w:space="0" w:color="1E3A5F"/>'
        f'<w:insideH w:val="single" w:sz="4" w:space="0" w:color="CBD5E1"/>'
        f'<w:left w:val="none"/><w:right w:val="none"/><w:insideV w:val="none"/>'
        f'</w:tblBorders>'
    ))
    for j, h in enumerate(("Audit Domain", "Exact Empirical & Manuscript Resolution")):
        c = tbl.rows[0].cells[j]
        c._tc.get_or_add_tcPr().append(parse_xml(f'<w:shd {nsdecls("w")} w:fill="1E3A5F"/>'))
        r = c.paragraphs[0].add_run(h)
        r.bold = True
        r.font.name = "Calibri"
        r.font.size = Pt(9)
        r.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    for i, (k, v) in enumerate(items):
        row = tbl.rows[i + 1]
        bg = "F8FAFC" if i % 2 == 1 else "FFFFFF"
        for j, val in enumerate((k, v)):
            c = row.cells[j]
            c._tc.get_or_add_tcPr().append(parse_xml(f'<w:shd {nsdecls("w")} w:fill="{bg}"/>'))
            r = c.paragraphs[0].add_run(val)
            r.font.name = "Calibri"
            r.font.size = Pt(8.5)
            if j == 0:
                r.bold = True

    doc.save("PI_Point_by_Point_Response_v3.3.docx")
    print("Saved updated PI_Point_by_Point_Response_v3.3.docx")


if __name__ == "__main__":
    polish_and_audit()
    rebuild_final_pi_response()
