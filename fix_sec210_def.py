import re
import docx
from docx.shared import Pt

DOCX_PATH = r"D:\noncoding mutation finding surver\NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx"
doc = docx.Document(DOCX_PATH)

updated = False
for p in doc.paragraphs:
    if "For continuous threshold-free evaluation (ROC-AUC and PR-AUC), NC-MutScan v3.3 defines" in p.text:
        # Detect the exact polyA threshold stated earlier in Section 2.10 (if numeric)
        polya_m = re.search(r"poly(?:adenylation|A)[^.]{0,60}?([0-9]+\.[0-9]+)\s*bits", p.text, flags=re.IGNORECASE)
        polya_clause = (
            f"τ_polyA = {polya_m.group(1)} bits for canonical AAUAAA/AUUAAA polyadenylation signal loss"
            if polya_m
            else "canonical AAUAAA/AUUAAA polyadenylation signal loss (|ΔS_polyA|)"
        )

        new_def = (
            "For continuous threshold-free evaluation (ROC-AUC and PR-AUC), NC-MutScan v3.3 defines a unified rule-integrated score "
            "S_unified (in bits) across all active regulatory modules after multi-TF competition and AID/SHM WRCY/RGYW passenger masking: "
            f"each module is first evaluated against its own calling threshold τ_m (τ_TF = 4.2 bits for transcription-factor motif shift "
            f"|ΔPWM_TF|; τ_splice = 2.8 bits for U1/U2 splice donor/acceptor loss and 3.8 bits for cryptic splice gain; and {polya_clause}), "
            "and variants passing their respective module threshold τ_m are rescaled onto a common post-threshold scale as "
            "S_unified = max(|ΔS_m|, 4.2) bits (anchoring all active splice, polyadenylation, and TF calls at or above the 4.2-bit "
            "unified decision boundary while preserving exact disruption magnitude above 4.2 bits), whereas variants failing their "
            "respective module threshold τ_m or masked as competitive/AID passengers receive S_unified = 0 so that the binary driver "
            "call corresponds without exception to S_unified ≥ 4.2 bits."
        )

        new_txt = re.sub(
            r"For continuous threshold-free evaluation \(ROC-AUC and PR-AUC\), NC-MutScan v3\.3 defines.*?4\.2 bits\.",
            new_def,
            p.text,
            flags=re.DOTALL,
        )
        fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
        fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else Pt(10.5)
        for r in list(p.runs):
            p._p.remove(r._r)
        nr = p.add_run(new_txt)
        nr.font.name = fn
        nr.font.size = fs
        updated = True
        break

assert updated, "Target sentence in Section 2.10 not found!"
doc.save(DOCX_PATH)

chk = docx.Document(DOCX_PATH)
for p in chk.paragraphs:
    if "S_unified = max(|ΔS_m|, 4.2)" in p.text:
        print("[PASS] Verified Section 2.10 definition on disk:\n")
        print(p.text)
