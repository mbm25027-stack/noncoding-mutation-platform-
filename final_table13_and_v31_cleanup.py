import os
import re
import docx

TARGET = r"D:\noncoding mutation finding surver\NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx"
if not os.path.isfile(TARGET):
    TARGET = os.path.join(os.path.expanduser("~"), "Downloads", "NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx")

doc = docx.Document(TARGET)

def set_para_text(p, new_text, bold=None):
    fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
    fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else None
    was_bold = p.runs[0].bold if p.runs else False
    for r in list(p.runs):
        p._p.remove(r._r)
    nr = p.add_run(new_text)
    nr.font.name = fn
    if fs is not None:
        nr.font.size = fs
    nr.bold = was_bold if bold is None else bold
    nr.font.highlight_color = None

def copy_table_to_end(doc, src_tbl):
    new_tbl = doc.add_table(rows=len(src_tbl.rows), cols=len(src_tbl.columns))
    new_tbl.style = src_tbl.style
    for r_idx, row in enumerate(src_tbl.rows):
        for c_idx, cell in enumerate(row.cells):
            new_tbl.rows[r_idx].cells[c_idx].text = cell.text

# 1. Fix Section 3.9 heading and update Table 12/13 -> Supplementary Table S2/S3, Table 14 -> Table 12
t13_caption = None
for p in doc.paragraphs:
    txt = p.text.strip()
    if "v3.1 baseline suite" in txt:
        set_para_text(p, "3.9 Frozen NC-MutScan v3.3 algorithmic regression and latency verification", bold=True)
    elif txt.startswith("Table 13."):
        t13_caption = p
    elif not txt.startswith("Supplementary Table"):
        new_t = txt
        new_t = re.sub(r"(?<!Supplementary )\bTable 12\b", "Supplementary Table S2", new_t)
        new_t = re.sub(r"(?<!Supplementary )\bTable 13\b", "Supplementary Table S3", new_t)
        new_t = re.sub(r"\bTable 14\b", "Table 12", new_t)
        # Ensure no stray v3.1 remains outside historical prototype note in Section 2.10
        if "v3.1" in new_t and "prototypes (v3.0/v3.1)" not in new_t:
            new_t = new_t.replace("v3.1", "v3.3")
        if new_t != txt:
            set_para_text(p, new_t)

# 2. Find and physically relocate the 54-locus Table 13 to Supplementary Table S3 at the end
moved_t13 = False
for tbl in list(doc.tables):
    blob = "\n".join(c.text for r in tbl.rows for c in r.cells).lower()
    if ("intergenic" in blob and "wrcy" in blob) or ("54" in blob and "truncated" in blob):
        cap = doc.add_paragraph("Supplementary Table S3. Frozen NC-MutScan v3.3 54-locus GRCh38/hg38 runtime and canonical-motif verification summary (~39 µs/sequence).")
        set_para_text(cap, cap.text, bold=True)
        copy_table_to_end(doc, tbl)
        tbl._element.getparent().remove(tbl._element)
        if t13_caption is not None:
            t13_caption._element.getparent().remove(t13_caption._element)
        moved_t13 = True
        break

doc.save(TARGET)

# 3. Final complete audit
doc_final = docx.Document(TARGET)
print("Relocated 54-locus Table 13 -> Supplementary Table S3:", moved_t13)
print("Total tables remaining in document:", len(doc_final.tables))
v31_hits = [p.text[:100] for p in doc_final.paragraphs if "v3.1" in p.text]
print("Remaining 'v3.1' mentions (should only be historical Section 2.10):", len(v31_hits))
for h in v31_hits:
    print("  ->", h)
for p in doc_final.paragraphs:
    if re.match(r"^(Table \d+\.|Supplementary Table S\d+\.|3\.9 |5\. Conclusions)", p.text.strip()):
        print("Caption/Sec:", p.text.strip()[:130])
