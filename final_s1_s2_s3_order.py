import re
import docx

TARGET = r"D:\noncoding mutation finding surver\NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx"
doc = docx.Document(TARGET)
body = doc._body._body

# 1. Sync any legacy 2.00 ms / 37 µs / 39 µs runtime text inside Supplementary Table S3 cells
s2_tbl = None
s3_tbl = None
for tbl in doc.tables:
    blob = "\n".join(c.text for r in tbl.rows for c in r.cells)
    if "TC-01" in blob:
        s2_tbl = tbl
    elif ("54" in blob or "Intergenic" in blob or "Canonical" in blob) and ("ms" in blob or "µs" in blob or "us" in blob):
        s3_tbl = tbl
        for row in tbl.rows:
            for c in row.cells:
                c.text = re.sub(r"\b2\.00\s*ms\b", "26.77 ms", c.text)
                c.text = re.sub(r"\b3[79](?:\.\d+)?\s*[µu]s\b", "495.7 µs", c.text)
                c.text = re.sub(r"~?2[57],\d{3}\s*seq(?:uences)?/s", "~2,017 sequences/s", c.text)

# 2. Locate Supplementary Table S1's footnote paragraph as the exact anchor
anchor_el = None
for p in doc.paragraphs:
    if p.text.strip().startswith("Supplementary Table S1 reports exact Average Precision"):
        anchor_el = p._p
        break

# 3. Move Supplementary Note S2, Caption S2, Table S2, Caption S3, Table S3 immediately after anchor_el
if anchor_el is not None:
    note_p = cap_s2 = cap_s3 = None
    for p in doc.paragraphs:
        st = p.text.strip()
        if st.startswith("Supplementary Note S2."):
            note_p = p._p
        elif st.startswith("Supplementary Table S2."):
            cap_s2 = p._p
        elif st.startswith("Supplementary Table S3."):
            cap_s3 = p._p

    ordered_els = [el for el in (note_p, cap_s2, s2_tbl._element if s2_tbl else None, cap_s3, s3_tbl._element if s3_tbl else None) if el is not None]
    cur = anchor_el
    for el in ordered_els:
        body.remove(el)
        cur.addnext(el)
        cur = el

doc.save(TARGET)

d_chk = docx.Document(TARGET)
for p in d_chk.paragraphs:
    if p.text.strip().startswith(("Supplementary Note", "Supplementary Table")):
        print("Final order:", p.text.strip()[:120])
