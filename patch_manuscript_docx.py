"""
patch_manuscript_docx.py
Directly updates the manuscript .docx file in-place to resolve all 8 PI review concerns:
  1. Unifies pipeline/software version references to v3.3.
  2. Purges hg19 chr7:140,453,135-140,453,136 (BRAF p.V600E coding artifact) and harmonizes
     coordinates to GRCh38 (with hg19 cross-references).
  3. Replaces Table 2 with the full N=150 locked blind chr12-22 benchmark, architectural
     ablations (Splice/PolyA-only, TF-only, Unpaired PWM), and SpliceAI (DS_max >= 0.20) comparison.
  4. Replaces Table 4 with the independent PCAWG Skin-Melanoma WGS (N=70 donors, 7,987,562 SNVs)
     donor-level zero-included NB validation, reclassifying AL358472.4, AC040160.1, and TCF3
     as recurrent UV-associated core-promoter ETS-motif hotspots.
  5. Replaces Table 6 with the Histology-Adjusted (+ histology_abbreviation) and Within-B-Cell-
     Lymphoma PCAWG cis-eQTL OLS results (confirming BCL2 beta=+0.897, qBH=4.56e-15 while
     reporting lineage attenuation of BCL6 qBH=0.210 and PAX5 qBH=0.758).
  6. Replaces "PPV" on 1:1 balanced cohorts with "Balanced Precision (50% prevalence)" plus
     Bayes-adjusted PPVs at 1%, 5%, and 10% clinical priors.
  7. Inserts/updates the revised Results (cis-eQTL lineage sensitivity) and Discussion (WXS off-target
     scope, UV-damaged ETS hotspots, and triage positioning) subsections.
"""

import copy
import glob
import os
import re
import shutil
import subprocess
import sys

try:
    import docx
    from docx.shared import Pt, RGBColor, Inches
    from docx.oxml import OxmlElement, parse_xml
    from docx.oxml.ns import nsdecls, qn
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "python-docx"])
    import docx
    from docx.shared import Pt, RGBColor, Inches
    from docx.oxml import OxmlElement, parse_xml
    from docx.oxml.ns import nsdecls, qn


def find_manuscript_docx() -> list[str]:
    candidates = []
    for root, _, files in os.walk("."):
        if any(skip in root.lower() for skip in ("venv", ".git", "__pycache__", "node_modules")):
            continue
        for fn in files:
            if fn.lower().endswith(".docx") and not fn.startswith("~$") and "_v3.3_PI_Revised" not in fn:
                candidates.append(os.path.join(root, fn))
    return sorted(candidates)


def set_cell_shading(cell, hex_color: str):
    shading = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{hex_color}"/>')
    cell._tc.get_or_add_tcPr().append(shading)


def set_cell_margins(cell, top=60, bottom=60, left=90, right=90):
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


def build_styled_table(doc, headers: list[str], rows: list[list[str]]):
    tbl = doc.add_table(rows=1 + len(rows), cols=len(headers))
    tbl.autofit = True
    set_table_borders(tbl)

    # Header row
    hdr_cells = tbl.rows[0].cells
    trPr = tbl.rows[0]._tr.get_or_add_trPr()
    trPr.append(OxmlElement("w:tblHeader"))
    trPr.append(OxmlElement("w:cantSplit"))
    for j, h_text in enumerate(headers):
        c = hdr_cells[j]
        set_cell_shading(c, "1E3A5F")
        set_cell_margins(c, top=80, bottom=80, left=80, right=80)
        p = c.paragraphs[0]
        p.paragraph_format.space_before = Pt(0)
        p.paragraph_format.space_after = Pt(0)
        run = p.add_run(h_text)
        run.bold = True
        run.font.name = "Calibri"
        run.font.size = Pt(8.5)
        run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    # Data rows
    for i, row_vals in enumerate(rows):
        row = tbl.rows[i + 1]
        r_trPr = row._tr.get_or_add_trPr()
        r_trPr.append(OxmlElement("w:cantSplit"))
        bg = "F8FAFC" if i % 2 == 1 else "FFFFFF"
        is_highlight = any(k in row_vals[0] for k in ("NC-MutScan v3.3", "TERT", "BCL2"))
        for j, val in enumerate(row_vals):
            c = row.cells[j]
            set_cell_shading(c, "EFF6FF" if is_highlight else bg)
            set_cell_margins(c, top=55, bottom=55, left=80, right=80)
            p = c.paragraphs[0]
            p.paragraph_format.space_before = Pt(0)
            p.paragraph_format.space_after = Pt(0)
            run = p.add_run(val)
            run.font.name = "Calibri"
            run.font.size = Pt(8.5)
            run.font.color.rgb = RGBColor(0x0F, 0x17, 0x2A)
            if j == 0 or is_highlight and j in (0, 4, 6, 7):
                run.bold = True
    return tbl


def replace_paragraph_text_preserve_style(p, new_text: str):
    font_name = "Calibri"
    font_size = Pt(11)
    is_bold = False
    is_italic = False
    if p.runs:
        r0 = p.runs[0]
        if r0.font.name:
            font_name = r0.font.name
        if r0.font.size:
            font_size = r0.font.size
        is_bold = bool(r0.bold)
        is_italic = bool(r0.italic)
    for r in list(p.runs):
        p._p.remove(r._r)
    nr = p.add_run(new_text)
    nr.font.name = font_name
    nr.font.size = font_size
    nr.bold = is_bold
    nr.italic = is_italic


def patch_docx(path: str) -> str:
    doc = docx.Document(path)
    stats = {
        "paragraphs_updated": 0,
        "braf_artifacts_removed": 0,
        "tables_replaced": 0,
        "sections_updated": 0,
    }

    # Regex substitutions applied across prose paragraphs
    text_rules = [
        (re.compile(r"\bNC-MutScan\s+v3\.[012]\b", re.I), "NC-MutScan v3.3"),
        (re.compile(r"\bNC-DriverEngine\s+v3\.[012]\b", re.I), "NC-DriverEngine v3.3"),
        (re.compile(r"\b97\.37%\s*PPV\b", re.I), "97.37% Balanced Precision (50% prevalence; Bayes-adjusted PPV: 27.3%, 66.1%, and 80.4% at 1%, 5%, and 10% clinical priors)"),
        (re.compile(r"\bPPV\s*=\s*97\.37%\b", re.I), "Balanced Precision (50% prevalence) = 97.37% (Bayes-adjusted PPV: 27.3%, 66.1%, and 80.4% at 1%, 5%, and 10% clinical priors)"),
        (re.compile(r"\bPositive Predictive Value \(PPV\)\b", re.I), "Balanced Precision (50% prevalence)"),
        (re.compile(r"\bchr7:140,?453,?13[56]\b"), "GRCh38 chr7:140,753,336 (hg19 chr7:140,453,136; BRAF p.V600E coding exon, excluded)"),
        (re.compile(r"\b140,?453,?13[56]\b"), "GRCh38 chr7:140,753,336 (excluded BRAF p.V600E coding locus)"),
    ]

    # 1. Patch prose paragraphs
    for p in list(doc.paragraphs):
        orig = p.text
        if not orig.strip():
            continue
        updated = orig
        for pat, repl in text_rules:
            updated = pat.sub(repl, updated)

        # Soften overclaimed driver phrasing for AC040160.1 / AL358472.4 / TCF3
        if any(g in updated for g in ("AC040160.1", "AL358472.4", "TCF3")) and "driver" in updated.lower():
            updated = re.sub(
                r"\b(novel|confirmed|putative)\s+noncoding\s+driver(s)?\b",
                r"recurrent UV-associated core-promoter ETS-motif hotspot\2",
                updated,
                flags=re.I,
            )

        # Flag lineage attenuation if an old paragraph cites BCL6 q=1.68e-11 or PAX5 q=5.72e-3 as standalone drivers
        if ("BCL6" in updated and ("1.68" in updated or "4.20" in updated or "0.903" in updated)):
            updated += (
                " [Revision Note: After adding tumor histology fixed effects (+ histology_abbreviation) "
                "and restricting within B-cell lymphomas (Lymph-BNHL + Lymph-CLL), BCL6 (beta = +0.075, "
                "qBH = 0.210) and PAX5 (beta = +0.013, qBH = 0.758) attenuated to non-significance due to "
                "baseline B-cell lineage expression, whereas BCL2 remained highly significant (beta = +0.897, "
                "qBH = 4.56e-15; within Lymph-BNHL alone: beta = +1.000, qBH = 1.48e-10).]"
            )

        if updated != orig:
            replace_paragraph_text_preserve_style(p, updated)
            stats["paragraphs_updated"] += 1

    # 2. Define canonical replacement tables
    t2_headers = [
        "Model / Variant Scoring Configuration", "Scored (N)", "TP", "FP", "TN", "FN",
        "Sensitivity (95% CI)", "Specificity (95% CI)", "Balanced Precision (50% Prev.)", "F1 Score"
    ]
    t2_rows = [
        ["NC-MutScan v3.3 (Full Differential ΔPWM)", "150", "37", "1", "74", "38", "49.33% (38.3–60.4%)", "98.67% (92.8–99.8%)", "97.37%", "0.6549"],
        ["Ablation 1: Splice / PolyA-PWM Only (no TFs)", "150", "36", "0", "75", "39", "48.00% (37.1–59.1%)", "100.00% (95.1–100.0%)", "100.00%", "0.6486"],
        ["Ablation 2: TF-Motif-PWM Only (no Splice/PolyA)", "150", "1", "1", "74", "74", "1.33% (0.2–7.2%)", "98.67% (92.8–99.8%)", "50.00%", "0.0260"],
        ["Ablation 3: Unpaired Single-Sequence PWM", "150", "65", "66", "9", "10", "86.67% (77.2–92.6%)", "12.00% (6.4–21.3%)", "49.62%", "0.6311"],
        ["External Comparator: SpliceAI (DS_max ≥ 0.20, GRCh38)", "139", "63", "3", "67", "6", "91.30% (82.3–96.0%)", "95.71% (88.1–98.5%)", "95.45%", "0.9333"],
    ]

    t4_headers = [
        "Locus / Gene", "GRCh38 Coordinates (hg19 Cross-Ref)", "Mutated WGS Donors (N=70)",
        "Total SNVs", "UV C>T/G>A", "Donor FE", "Donor pNB", "Genome-Wide qNB", "Revised Classification"
    ]
    t4_rows = [
        ["TERT (C228T/C250T)", "chr5:1,295,113–1,295,135 (hg19 chr5:1,295,228–1,295,250)", "24 / 70 (34.3%)", "25", "19/25 (76.0%)", "178.8x", "4.44e-16", "2.53e-08", "Confirmed noncoding driver (ETS gain)"],
        ["AL358472.4 promoter", "chr1:153,990,702–153,990,745 (hg19 chr1:153,963,178–153,963,221)", "13 / 70 (18.6%)", "15", "14/15 (93.3%)", "96.9x", "1.28e-11", "7.27e-04", "Recurrent UV-associated ETS hotspot"],
        ["AC040160.1 promoter", "chr16:67,227,004–67,227,077 (hg19 chr16:67,260,907–67,260,980)", "12 / 70 (17.1%)", "17", "14/17 (82.4%)", "89.4x", "8.39e-11", "4.78e-03", "Recurrent UV-associated ETS hotspot"],
        ["TCF3 promoter", "chr19:1,652,862–1,653,023 (hg19 chr19:1,652,862–1,653,023)", "11 / 70 (15.7%)", "15", "15/15 (100.0%)", "82.0x", "5.52e-10", "3.14e-02", "Recurrent UV-associated ETS hotspot"],
    ]

    t6_headers = [
        "Target Gene", "Locus (GRCh38 / hg19)", "OLS Specification", "WT / MUT",
        "Beta_mut", "95% CI", "Two-Sided P", "BH qBH", "Biological Interpretation"
    ]
    t6_rows = [
        ["BCL2", "chr18:63,318,567 (hg19:60,985,800)", "Within B-Cell Lymphomas (Lymph-BNHL+CLL)", "120 / 46", "+0.897", "[+0.698, +1.096]", "1.14e-15", "4.56e-15", "True cis-regulatory driver (lineage-independent)"],
        ["BCL2", "chr18:63,318,567 (hg19:60,985,800)", "Within Lymph-BNHL Alone (GC B-NHL)", "55 / 43", "+1.000", "[+0.735, +1.264]", "3.71e-11", "1.48e-10", "True cis-regulatory driver (within GC B-NHL)"],
        ["BCL6", "chr3:187,745,212 (hg19:187,463,000)", "Histology-Adjusted (+ histology fixed effects)", "263 / 96", "+0.075", "[-0.016, +0.165]", "1.05e-01", "2.10e-01", "Attenuates after lineage adjustment (AID SHM target)"],
        ["BCL6", "chr3:187,745,212 (hg19:187,463,000)", "Within B-Cell Lymphomas (Lymph-BNHL+CLL)", "75 / 91", "+0.051", "[-0.062, +0.163]", "3.75e-01", "7.43e-01", "Non-significant within matched B-cell lymphomas"],
        ["PAX5", "chr9:37,034,003 (hg19:37,034,000)", "Within B-Cell Lymphomas (Lymph-BNHL+CLL)", "142 / 24", "+0.013", "[-0.065, +0.090]", "7.43e-01", "7.58e-01", "Attenuates after lineage adjustment (AID SHM target)"],
        ["WDR74", "chr11:62,841,608 (hg19:62,609,080)", "Histology-Adjusted (+ histology fixed effects)", "398 / 27", "-0.009", "[-0.070, +0.051]", "7.58e-01", "7.58e-01", "Confirmed passenger / negative control"],
    ]

    # 3. Inspect and upgrade tables in-place
    replaced_t2, replaced_t4, replaced_t6 = False, False, False
    for tbl in list(doc.tables):
        flat_text = " ".join(c.text.strip() for row in tbl.rows for c in row.cells)

        # Purge any BRAF V600E row in any remaining table
        for row in list(tbl.rows)[1:]:
            row_txt = " ".join(c.text.strip() for c in row.cells)
            if "140453135" in row_txt or "140,453,135" in row_txt or "140453136" in row_txt or "140,453,136" in row_txt:
                tbl._tbl.remove(row._tr)
                stats["braf_artifacts_removed"] += 1

        # Match Blind Benchmark Table (Table 2)
        if not replaced_t2 and ("98.67" in flat_text or "49.33" in flat_text or ("Sensitivity" in flat_text and "Specificity" in flat_text)):
            new_tbl = build_styled_table(doc, t2_headers, t2_rows)
            tbl._element.addprevious(new_tbl._element)
            tbl._element.getparent().remove(tbl._element)
            replaced_t2 = True
            stats["tables_replaced"] += 1
            continue

        # Match Melanoma Core-Promoter Table (Table 4)
        if not replaced_t4 and ("AC040160" in flat_text or "AL358472" in flat_text or ("TERT" in flat_text and "TCF3" in flat_text)):
            new_tbl = build_styled_table(doc, t4_headers, t4_rows)
            tbl._element.addprevious(new_tbl._element)
            tbl._element.getparent().remove(tbl._element)
            replaced_t4 = True
            stats["tables_replaced"] += 1
            continue

        # Match PCAWG eQTL Table (Table 6)
        if not replaced_t6 and ("BCL6" in flat_text and "BCL2" in flat_text and ("PAX5" in flat_text or "WDR74" in flat_text)):
            new_tbl = build_styled_table(doc, t6_headers, t6_rows)
            tbl._element.addprevious(new_tbl._element)
            tbl._element.getparent().remove(tbl._element)
            replaced_t6 = True
            stats["tables_replaced"] += 1
            continue

    # 4. Append/Insert Revised Sensitivity & Scope Subsections before References (or at end)
    ref_para = None
    for p in doc.paragraphs:
        if p.text.strip().lower() in ("references", "bibliography", "6. references", "7. references"):
            ref_para = p
            break

    def add_section_block(anchor_p, heading_text: str, body_text: str):
        if anchor_p is not None:
            hp = anchor_p.insert_paragraph_before(heading_text)
            bp = anchor_p.insert_paragraph_before(body_text)
        else:
            hp = doc.add_paragraph(heading_text)
            bp = doc.add_paragraph(body_text)
        for r in hp.runs:
            r.bold = True
            r.font.name = "Calibri"
            r.font.size = Pt(12)
            r.font.color.rgb = RGBColor(0x1E, 0x3A, 0x5F)
        for r in bp.runs:
            r.font.name = "Calibri"
            r.font.size = Pt(11)
        bp.paragraph_format.space_after = Pt(8)
        stats["sections_updated"] += 1
        return bp

    eqtl_body = (
        "To ensure that multi-cohort cis-eQTL associations were not confounded by baseline tissue-lineage expression "
        "differences, we augmented our multivariable ordinary least squares (OLS) model "
        "(log2[FPKM-UQ + 1] ~ I_mut + CNV + Purity + log10[TMB + 1]) with tumor histology fixed effects "
        "(+ histology_abbreviation) and performed stratified regressions restricted strictly within B-cell lymphoid "
        "malignancies (Lymph-BNHL and Lymph-CLL, N = 166 matched donors). Under both histology-adjusted and "
        "within-lymphoma specifications, promoter/5'-UTR mutations in BCL2 (WT = 120, MUT = 46) remained strongly "
        "associated with elevated BCL2 transcript abundance (beta = +0.897, 95% CI [+0.698, +1.096], p = 1.14e-15, "
        "qBH = 4.56e-15), persisting when restricted further to germinal-center B-cell non-Hodgkin lymphoma alone "
        "(Lymph-BNHL, WT = 55, MUT = 43, beta = +1.000, 95% CI [+0.735, +1.264], qBH = 1.48e-10). Conversely, "
        "apparent pan-cohort expression elevations at BCL6 (beta = +0.075, qBH = 0.210 after histology adjustment; "
        "beta = +0.051, qBH = 0.743 within B-cell lymphomas) and PAX5 (beta = +0.013, qBH = 0.758 within B-cell "
        "lymphomas) attenuated to statistical non-significance once baseline B-cell lineage expression was separated "
        "from non-lymphoid carcinomas, consistent with activation-induced cytidine deaminase (AID) somatic "
        "hypermutation targeting actively transcribed germinal-center lineage genes rather than inducing steady-state "
        "transcriptional gain. Consistent with prior PCAWG passenger benchmarks, WDR74 promoter mutations showed no "
        "expression effect (beta = -0.009, qBH = 0.758)."
    )

    disc_body = (
        "Several methodological boundaries define the scope of NC-MutScan v3.3. First, in whole-exome cohorts such as "
        "TCGA-SKCM (N = 466 donors with somatic MAF records out of 470 clinical cases), core-promoter detection relies "
        "on off-target flanking reads extending upstream of exon-1 capture baits; while this enables opportunistic "
        "promoter screening in legacy WXS datasets, fold-enrichment estimates relative to exome-wide mutation rates "
        "are sensitive to capture depth gradients. Independent validation in unbiased whole-genome sequencing "
        "(PCAWG Skin-Melanoma / MELA-AU, N = 70 donors, 7,987,562 SNVs) confirmed high-frequency donor recurrence "
        "across TERT (24/70, 34.3%), AL358472.4 (13/70, 18.6%), AC040160.1 (12/70, 17.1%), and TCF3 (11/70, 15.7%). "
        "Because 82.4%–100.0% of SNVs at AL358472.4, AC040160.1, and TCF3 represent C>T/G>A pyrimidine-dimer "
        "transitions within CTTCCGG ETS-binding motifs—a signature of impaired nucleotide excision repair at "
        "ETS-occupied promoters—we classify these three loci as recurrent UV-associated core-promoter ETS-motif "
        "hotspots pending matched transcriptomic or reporter validation rather than confirmed oncogenic drivers. "
        "Second, on the locked blind chr12–22 benchmark (N = 150; calibrated on N = 955 chr1–11, X variants), paired "
        "wild-type versus mutant ΔPWM scoring improved specificity by +86.7 percentage points over unpaired "
        "single-sequence motif scanning (98.67% vs. 12.00%) and achieved a >3-fold lower false-positive rate "
        "(1.33% vs. 4.29%) than SpliceAI (DS_max >= 0.20, N = 139 scored), albeit with lower recall (49.33% vs. "
        "91.30%) on deep-intronic splice variants lacking canonical PWM consensus matches. NC-MutScan v3.3 is "
        "therefore best positioned as a zero-dependency, ultra-high-specificity triage filter for rapid somatic and "
        "clinical variant prioritization."
    )

    add_section_block(
        ref_para,
        "Histology-Adjusted and Lineage-Stratified cis-eQTL Sensitivity Analysis",
        eqtl_body,
    )
    if not replaced_t6:
        t6_new = build_styled_table(doc, t6_headers, t6_rows)
        if ref_para is not None:
            ref_para._element.addprevious(t6_new._element)

    add_section_block(
        ref_para,
        "Methodological Scope, Off-Target WXS Discovery, and UV-Associated ETS Hotspots",
        disc_body,
    )
    if not replaced_t2:
        t2_new = build_styled_table(doc, t2_headers, t2_rows)
        if ref_para is not None:
            ref_para._element.addprevious(t2_new._element)
    if not replaced_t4:
        t4_new = build_styled_table(doc, t4_headers, t4_rows)
        if ref_para is not None:
            ref_para._element.addprevious(t4_new._element)

    base, ext = os.path.splitext(path)
    out_path = f"{base}_v3.3_PI_Revised{ext}"
    doc.save(out_path)
    print(f"Patched: {path}")
    print(f"Saved  : {out_path}")
    print(f"Summary: {stats} | Tables replaced in-place: T2={replaced_t2}, T4={replaced_t4}, T6={replaced_t6}")
    return out_path


def main():
    files = find_manuscript_docx()
    if not files:
        print("No .docx files found in workspace! Listing current directory:")
        print(os.listdir("."))
        return
    print("Found .docx files:", files)
    for f in files:
        patch_docx(f)


if __name__ == "__main__":
    main()
