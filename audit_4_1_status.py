import os
import re
import docx

TARGET = r"D:\noncoding mutation finding surver\NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx"
if not os.path.isfile(TARGET):
    TARGET = os.path.join(os.path.expanduser("~"), "Downloads", "NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx")

doc = docx.Document(TARGET)
full_txt = "\n".join(p.text for p in doc.paragraphs)

has_zero_auc = "ROC-AUC = 0.000" in full_txt or "ROC=0.000" in full_txt
print("Target checked:", TARGET)
print("Contains stale 0.000 AUC values?:", has_zero_auc)

for p in doc.paragraphs:
    if "ROC-AUC =" in p.text or "Figure 18." in p.text:
        print("  ->", p.text[:180])

# Scrub any residual 100%/100% or 54/54 accuracy headlines in main text
changed = 0
for p in doc.paragraphs:
    t = p.text
    if any(k in t for k in ("100% sensitivity and 100% specificity", "54/54", "50/50")):
        if not t.startswith("Supplementary Note S2") and not t.startswith("Table 12") and not t.startswith("Table 13"):
            new_t = (
                t.replace("100% sensitivity and 100% specificity", "deterministic rule concordance (Supplementary Note S2)")
                 .replace("54/54 (100%)", "all 54 verification loci (Supplementary Note S2)")
                 .replace("54/54", "all 54 verification loci (Supplementary Note S2)")
                 .replace("50/50 (100%)", "all 50 canonical-motif unit tests (Supplementary Note S2)")
                 .replace("50/50", "all 50 canonical-motif unit tests (Supplementary Note S2)")
            )
            if new_t != t:
                fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
                fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else None
                for r in list(p.runs):
                    p._p.remove(r._r)
                nr = p.add_run(new_t)
                nr.font.name = fn
                if fs is not None:
                    nr.font.size = fs
                changed += 1

doc.save(TARGET)
print(f"Scrubbed {changed} residual self-referential headline lines in {os.path.basename(TARGET)}.")
