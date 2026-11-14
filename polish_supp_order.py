import docx

fp = r"D:\noncoding mutation finding surver\NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx"
doc = docx.Document(fp)

def set_p(p, txt, bold=True):
    fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
    fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else None
    for r in list(p.runs):
        p._p.remove(r._r)
    nr = p.add_run(txt)
    nr.font.name = fn
    if fs:
        nr.font.size = fs
    nr.bold = bold

# 1. Align Table 5 caption with TERT positive-control + UV-hotspot framing
for p in doc.paragraphs:
    if p.text.strip().startswith("Table 5. Prioritized core-promoter driver candidates"):
        set_p(p, "Table 5. Prioritized core-promoter recurrent hotspots (≤300 bp from the TSS, GRCh38/hg38; TERT literature positive control and UV-associated ETS hotspots)")

# 2. Move Supplementary Table S2 & S3 blocks after Supplementary Table S1 so S1 -> S2 -> S3 is sequential
body = doc._body._body
s1_tbl_el = None
for tbl in doc.tables:
    if "LOOCV" in "\n".join(c.text for r in tbl.rows for c in r.cells):
        s1_tbl_el = tbl._element
        break

if s1_tbl_el is not None:
    move_els = []
    for p in doc.paragraphs:
        if p.text.strip().startswith(("Supplementary Note S2.", "Supplementary Table S2.", "Supplementary Table S3.")):
            move_els.append(p._p)
    for tbl in doc.tables:
        blob = "\n".join(c.text for r in tbl.rows for c in r.cells)
        if "TC-01" in blob or "Intergenic background" in blob:
            move_els.append(tbl._element)
    anchor = s1_tbl_el
    for el in move_els:
        body.remove(el)
        anchor.addnext(el)
        anchor = el

doc.save(fp)
print("Aligned Table 5 caption and ordered Supplementary Tables S1 -> S2 -> S3.")
