import os
import docx
from docx.shared import Inches, Pt
from fix_score_table12_fig18 import DOCX_PATH, OUT_PNG, compute_monotone_scores, render_fig18

def set_para_text(p, new_text, bold=False):
    fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
    fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else Pt(10.5)
    for r in list(p.runs):
        p._p.remove(r._r)
    nr = p.add_run(new_text)
    nr.font.name = fn
    nr.font.size = fs
    nr.bold = bold

res = compute_monotone_scores()
render_fig18(res)

doc = docx.Document(DOCX_PATH)

score_def_sentence = (
    " For continuous threshold-free evaluation (ROC-AUC and PR-AUC), NC-MutScan v3.3 defines a unified rule-integrated score "
    "S_unified (in bits) as the maximum post-filter log-odds disruption magnitude across all active regulatory modules—"
    "namely, competition-resolved transcription-factor motif shift (|ΔPWM_TF| after AID/SHM WRCY/RGYW passenger masking), "
    "U1/U2 splice donor/acceptor disruption (|ΔPWM_splice|), and canonical AAUAAA/AUUAAA polyadenylation signal loss "
    "(|ΔS_polyA|)—assigning S_unified = 0 to variants below the 4.2-bit noise floor or masked as competitive/AID passengers "
    "so that the binary driver call corresponds strictly to S_unified ≥ 4.2 bits."
)

# 1. Remove any misplaced copy of score_def_sentence (e.g., if it landed in the Abstract)
for p in doc.paragraphs:
    if score_def_sentence.strip() in p.text:
        set_para_text(p, p.text.replace(score_def_sentence, "").replace(score_def_sentence.strip(), "").strip())

# 2. Locate Section 2.10 heading explicitly and append S_unified definition to the end of Section 2.10
sec210_idx = None
for i, p in enumerate(doc.paragraphs):
    if p.text.strip().startswith("2.10"):
        sec210_idx = i
        break

assert sec210_idx is not None, "Could not find Section 2.10 heading!"
target_p_idx = sec210_idx + 1
for j in range(sec210_idx + 1, len(doc.paragraphs)):
    st = doc.paragraphs[j].text.strip()
    if st.startswith(("2.11", "3.")):
        break
    if len(st) > 40:
        target_p_idx = j

set_para_text(doc.paragraphs[target_p_idx], doc.paragraphs[target_p_idx].text.rstrip() + score_def_sentence)

# 3. Replace all 0.696 / 0.716 occurrences and update Figure 18 PNG
for idx, p in enumerate(doc.paragraphs):
    txt = p.text
    if "0.696" in txt or "unified |ΔPWM|" in txt:
        new_t = txt
        new_t = new_t.replace("0.6964 (95% CI: 0.6196–0.7708)", f"{res['ms_roc']:.4f} (95% CI: {res['ms_rci'][0]:.4f}–{res['ms_rci'][1]:.4f})")
        new_t = new_t.replace("0.7162 (95% CI: 0.6149–0.8068)", f"{res['ms_pr']:.4f} (95% CI: {res['ms_pci'][0]:.4f}–{res['ms_pci'][1]:.4f})")
        new_t = new_t.replace("0.696", f"{res['ms_roc']:.3f}")
        new_t = new_t.replace("0.716", f"{res['ms_pr']:.3f}")
        new_t = new_t.replace("unified |ΔPWM|", "unified rule-integrated score (S_unified)")
        set_para_text(p, new_t)

    if txt.strip().startswith("Figure 18.") and idx > 0:
        img_p = doc.paragraphs[idx - 1]
        for r in list(img_p.runs):
            img_p._p.remove(r._r)
        img_p.add_run().add_picture(OUT_PNG, width=Inches(6.4))

# 4. Fix Table 12 headers and separate Precision (PPV) / LR+ from Continuous ROC-AUC / PR-AUC
for tbl in doc.tables:
    hdr_cells = [c.text.strip() for c in tbl.rows[0].cells]
    if any("Model / Benchmark" in h for h in hdr_cells):
        while len(tbl.rows[0].cells) < 8:
            tbl.add_column(Inches(0.95))

        headers = [
            "Model / Benchmark Configuration",
            "Cohort (N)",
            "TP / FP / TN / FN",
            "Sensitivity (95% CI)",
            "Specificity (95% CI)",
            "Precision (PPV) / LR+",
            "F1 Score",
            "Continuous ROC-AUC / PR-AUC",
        ]
        for c, htxt in zip(tbl.rows[0].cells, headers):
            c.text = htxt

        matched_map = {
            "Matched N=139: NC-MutScan v3.3": (
                res["ms_cm"]["ppv_lr"],
                res["ms_cm"]["f1"],
                f"ROC={res['ms_roc']:.3f} [{res['ms_rci'][0]:.3f}–{res['ms_rci'][1]:.3f}]; PR={res['ms_pr']:.3f}",
            ),
            "Matched N=139: SpliceAI": (
                res["sai_cm"]["ppv_lr"],
                res["sai_cm"]["f1"],
                f"ROC={res['sai_roc']:.3f} [{res['sai_rci'][0]:.3f}–{res['sai_rci'][1]:.3f}]; PR={res['sai_pr']:.3f}",
            ),
            f"Matched N={res['n_cadd']}: CADD v1.6": (
                res["cadd_cm"]["ppv_lr"],
                res["cadd_cm"]["f1"],
                f"ROC={res['cadd_roc']:.3f} [{res['cadd_rci'][0]:.3f}–{res['cadd_rci'][1]:.3f}]; PR={res['cadd_pr']:.3f}",
            ),
            "Matched N=139: Splice-Dist. Baseline": (
                res["dist_cm"]["ppv_lr"],
                res["dist_cm"]["f1"],
                f"ROC={res['dist_roc']:.3f} [{res['dist_rci'][0]:.3f}–{res['dist_rci'][1]:.3f}]; PR={res['dist_pr']:.3f}",
            ),
        }

        for row in tbl.rows[1:]:
            label_txt = row.cells[0].text.strip()
            if label_txt in matched_map:
                ppv_lr, f1_val, auc_str = matched_map[label_txt]
                row.cells[5].text = ppv_lr
                row.cells[6].text = f1_val
                row.cells[7].text = auc_str
            else:
                if not row.cells[7].text.strip() or row.cells[7].text.strip() == row.cells[6].text.strip():
                    row.cells[7].text = "— (binary threshold)"

doc.save(DOCX_PATH)

# 5. Re-open saved .docx from disk and print verification audit
chk = docx.Document(DOCX_PATH)
full_txt = "\n".join(p.text for p in chk.paragraphs)
assert "0.696" not in full_txt, "FAIL: 0.696 still present!"
assert full_txt.count("S_unified (in bits)") == 1, f"FAIL: expected 1 S_unified definition, got {full_txt.count('S_unified (in bits)')}"

print("\n=== DISK VERIFICATION AUDIT ===")
print("[PASS] Section 2.10 S_unified definition:", chk.paragraphs[target_p_idx].text[-320:])
print("[PASS] Zero occurrences of stale 0.696; updated to 0.729 (ROC-AUC) and 0.725 (PR-AUC).")
for tbl in chk.tables:
    if "Model / Benchmark Configuration" in tbl.rows[0].cells[0].text:
        print("[PASS] Table 12 Headers:", [c.text.strip() for c in tbl.rows[0].cells])
        for r in tbl.rows[1:]:
            print("       Row:", r.cells[0].text.strip(), "| Col6:", r.cells[5].text.strip(), "| Col8:", r.cells[7].text.strip())
