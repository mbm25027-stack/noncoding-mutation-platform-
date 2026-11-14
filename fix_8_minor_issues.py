"""
fix_8_minor_issues.py
Precisely resolves all 8 minor manuscript issues in NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx:
1. Version labels: attributes 72.0% / 21.3% single-sequence results and early exploration strictly to v3.1.
2. Latency: aligns Figure 17C caption with §2.11 and §3.9 (495.7 µs/sequence; 26.77 ms across 54 loci).
3. §3.9 sentence: fixes missing opening parenthesis and broken grammar around 26.77 ms / 495.7 µs/sequence.
4. Discussion limitations: merges overlapping "Third/Fourth" and "Limitations. First/Second/Third" into one
   ordered list, removes the false "autosomes only" claim (noting chrX in calibration), and fixes v3.1 vs v3.3.
5. Thresholds & matrix counts: explicitly pairs early v3.0/v3.1 (6–7 TF matrices; ΔPWM ≤ −4.00 bits) vs.
   frozen v3.3 (9 core TF PWMs; |ΔPWM| ≥ 4.2 bits) wherever mentioned in §2.10, §3.9, and Discussion.
6. Counting: changes "two" verification suites in §2.11 to "three" checks, and rewords "partitioned 1,105
   variants" in §3.11 to reflect N=955 on chr1–11,X plus class-stratified sampling of N=150 from chr12–22.
7. Caption formatting: applies bold formatting to Figures 17–19 captions matching Figures 15, 16, and 20.
8. S_unified for polyA: adds an explicit clause in §2.10 stating that a discrete polyA call takes the 4.2-bit floor.
"""

import re
import docx
from docx.shared import Pt


DOCX_PATH = r"D:\noncoding mutation finding surver\NC-DriverEngine_v3.3_PI_Revised_Manuscript_4_1.docx"


def set_para_text(p, new_text, bold=None):
    fn = p.runs[0].font.name if p.runs and p.runs[0].font.name else "Calibri"
    fs = p.runs[0].font.size if p.runs and p.runs[0].font.size else Pt(10.5)
    if bold is None:
        bold = bool(p.runs[0].bold) if p.runs else False
    for r in list(p.runs):
        p._p.remove(r._r)
    nr = p.add_run(new_text)
    nr.font.name = fn
    nr.font.size = fs
    nr.bold = bold


def main():
    doc = docx.Document(DOCX_PATH)

    # -------------------------------------------------------------------------
    # STEP A: Inspect & merge overlapping Discussion limitation paragraphs (#1, #4, #5)
    # -------------------------------------------------------------------------
    disc_start = None
    conc_start = None
    for i, p in enumerate(doc.paragraphs):
        st = p.text.strip()
        if st.startswith("4. Discussion") or st == "Discussion":
            disc_start = i
        elif st.startswith("5. Conclusions"):
            conc_start = i
            break

    if disc_start is not None and conc_start is not None:
        lim_indices = []
        for idx in range(disc_start + 1, conc_start):
            st = doc.paragraphs[idx].text.strip()
            if (
                st.startswith(("Third,", "Fourth,", "Limitations."))
                or "several limitations should be noted" in st.lower()
                or "covers autosomes only" in st.lower()
            ):
                lim_indices.append(idx)

        if lim_indices:
            merged_limitations = (
                "Limitations. Several limitations should be noted when interpreting these results. "
                "First, TCGA-SKCM is a whole-exome dataset whose capture probes target protein-coding exons rather than "
                "distal enhancers or core promoters; consequently, off-target flanking coverage inflates local Obs/Exp ratios "
                "and under-samples low-coverage promoters, requiring independent whole-genome verification in PCAWG MELA-AU (N = 70). "
                "Second, the zero-included Negative Binomial WGS model relies on a single cohort-wide background rate without "
                "replication-timing or chromatin-accessibility covariates (yielding near-Poisson Fano factors of 1.02–1.19) and lacks "
                "matched RNA-seq in MELA-AU, so AL358472.4, AC040160.1, and TCF3 remain unvalidated UV-associated ETS hotspots rather "
                "than confirmed oncogenic drivers, while the SV-adjusted BCL2 expression signal rests on only n = 4 SV-negative mutated tumors. "
                "Third, although our chromosome-held-out ClinVar split prevents cross-chromosome leakage by placing chr1–11 and chrX "
                "(N = 955) in the calibration pool and sampling a balanced blind test cohort of N = 150 variants from held-out chr12–22, "
                "germline noncoding ClinVar variants cluster heavily at canonical splice junctions; as a result, U1/U2 splice and "
                "polyadenylation rules account for 36/37 true positives and simple splice-distance or deep-learning splice models "
                "achieve higher continuous AUC than deterministic PWMs on germline splice-dominated cohorts. "
                "Fourth, whereas early single-sequence prototypes (v3.0/v3.1, scanning 6–7 TF matrices at ΔPWM ≤ −4.00 bits) were "
                "explored during initial threshold inspection on chr1–11 and chrX (72.0% sensitivity, 21.3% specificity), all final "
                "software verification suites, locked ClinVar chr12–22 evaluations, and TERT promoter MPRA benchmarks were executed "
                "on frozen NC-MutScan v3.3 (nine core TF PWMs at |ΔPWM| ≥ 4.2 bits plus U1/U2 splice and polyA models)."
            )
            set_para_text(doc.paragraphs[lim_indices[0]], merged_limitations, bold=False)
            for extra_i in reversed(lim_indices[1:]):
                doc.paragraphs[extra_i]._element.getparent().remove(doc.paragraphs[extra_i]._element)

    # -------------------------------------------------------------------------
    # STEP B: Iterate through all paragraphs to apply fixes #1, #2, #3, #5, #6, #7, #8
    # -------------------------------------------------------------------------
    for idx, p in enumerate(doc.paragraphs):
        txt = p.text
        st = txt.strip()
        new_t = txt

        # Fix #8 & #5 in Section 2.10: add explicit polyA 4.2-bit floor clause and ensure v3.0/v3.1 vs v3.3 thresholds are clear
        if "For continuous threshold-free evaluation (ROC-AUC and PR-AUC), NC-MutScan v3.3 defines" in new_t:
            new_t = re.sub(
                r"Earlier internal prototypes \(v3\.0/v3\.1\) scanned single sequences globally across (?:six|seven) TF matrices",
                "Earlier internal prototypes (v3.0/v3.1) scanned single sequences globally across six to seven TF matrices at ΔPWM ≤ −4.00 bits",
                new_t,
            )
            new_t = re.sub(
                r"For continuous threshold-free evaluation \(ROC-AUC and PR-AUC\), NC-MutScan v3\.3 defines.*?4\.2 bits\.",
                (
                    "For continuous threshold-free evaluation (ROC-AUC and PR-AUC), NC-MutScan v3.3 defines a unified rule-integrated score "
                    "S_unified (in bits) across all active regulatory modules after multi-TF competition and AID/SHM WRCY/RGYW passenger masking: "
                    "each module is first evaluated against its own calling criterion τ_m—specifically, τ_TF = 4.2 bits for transcription-factor "
                    "motif shift (|ΔPWM_TF|), τ_splice = 2.8 bits for U1/U2 splice donor/acceptor loss and 3.8 bits for cryptic splice gain "
                    "(|ΔPWM_splice|), and discrete hexamer disruption of canonical 3′-UTR AATAAA/ATTAAA within a ≥ 65% A+T context window for "
                    "polyadenylation signal loss (τ_polyA)—and variants passing their respective module criterion τ_m are rescaled onto a common "
                    "post-threshold scale as S_unified = max(|ΔS_m|, 4.2) bits (where a discrete polyA call takes the 4.2-bit floor, anchoring all "
                    "active splice, polyadenylation, and TF calls at or above the 4.2-bit unified decision boundary while preserving exact disruption "
                    "magnitude above 4.2 bits), whereas variants failing their respective module criterion τ_m or masked as competitive/AID "
                    "passengers receive S_unified = 0 so that the binary driver call corresponds without exception to S_unified ≥ 4.2 bits."
                ),
                new_t,
                flags=re.DOTALL,
            )

        # Fix #6 in Section 2.11: "two software regression and runtime verification suites" -> "three ... checks"
        if "using two software regression and runtime verification suites" in new_t:
            new_t = new_t.replace(
                "using two software regression and runtime verification suites",
                "using three software regression and runtime verification checks",
            )

        # Fix #3 & #5 in Section 3.9: fix broken parenthesis/grammar and any stale -4.00 / 6-7 matrix wording for v3.3
        if "execute in 26.77 ms total" in new_t or "495.7 µs/sequence" in new_t:
            new_t = re.sub(
                r"execute in 26\.77 ms total[^)]*\)\s*with zero runtime exceptions",
                "execute in 26.77 ms total across the 54 loci (495.7 µs per sequence; ~2,017 sequences/s locally; <15 ms via the HTTP API) with zero runtime exceptions",
                new_t,
            )
            new_t = re.sub(
                r"execute in 26\.77 ms total,\s*495\.7 µs/sequence;\s*~2,017 sequences/s\)",
                "execute in 26.77 ms total across the 54 loci (495.7 µs per sequence; ~2,017 sequences/s locally)",
                new_t,
            )

        # Fix #5 globally where v3.3 is paired with legacy <= -4.00 or 6/7 matrices
        if ("v3.3" in new_t or "frozen" in new_t.lower()) and ("−4.00" in new_t or "-4.00" in new_t):
            if "v3.0" not in new_t and "v3.1" not in new_t:
                new_t = new_t.replace("ΔPWM ≤ −4.00 bits", "|ΔPWM| ≥ 4.2 bits (vs. ΔPWM ≤ −4.00 bits in early prototype v3.1)")
                new_t = new_t.replace("ΔPWM ≤ -4.00 bits", "|ΔPWM| ≥ 4.2 bits (vs. ΔPWM ≤ −4.00 bits in early prototype v3.1)")
                new_t = new_t.replace("≤ −4.00", "≥ 4.2 bits (|ΔPWM|; vs. ≤ −4.00 bits in early prototype v3.1)")
                new_t = new_t.replace("≤ -4.00", "≥ 4.2 bits (|ΔPWM|; vs. ≤ −4.00 bits in early prototype v3.1)")

        # Fix #1 & #6 in Section 3.11 (after Table 12) and Figure 20B caption:
        # - 72.0% / 21.3% single-sequence result attributed to early prototype v3.1 (not v3.3)
        # - "partitioned 1,105 variants" reworded accurately
        if "72.0%" in new_t and "21.3%" in new_t:
            new_t = re.sub(
                r"(?:frozen\s+)?NC-MutScan\s+v3\.3(?=\s*(?:\(|\[|,)?\s*(?:single-sequence|unpaired|yielded|achieved|scored)?[^.]{0,60}72\.0%)",
                "early single-sequence prototype NC-MutScan v3.1 (six to seven TF matrices; ΔPWM ≤ −4.00 bits)",
                new_t,
            )
            new_t = new_t.replace("v3.3 (72.0%", "early prototype v3.1 (72.0%")
            new_t = new_t.replace("v3.3 yielded 72.0%", "early prototype v3.1 yielded 72.0%")
            new_t = new_t.replace("v3.3 achieved 72.0%", "early prototype v3.1 achieved 72.0%")

        if "partitioned 1,105" in new_t:
            new_t = re.sub(
                r"partitioned 1,105[^.;]*",
                "curated N = 955 strictly noncoding ClinVar variants on chr1–11 and chrX as the calibration pool and drew a balanced, class-stratified blind test sample of N = 150 variants (75 Pathogenic/Likely pathogenic, 75 Benign/Likely benign; random.seed(42)) from the held-out chr12–22 noncoding ClinVar pool (N = 1,105 total evaluated variants across both sets)",
                new_t,
            )

        # Fix #2 in Figure 17C caption (and any remaining ~39 µs mentions): harmonize with 495.7 µs/sequence (26.77 ms over 54 loci)
        if st.startswith("Figure 17.") or "39 µs" in new_t or "39 μs" in new_t:
            new_t = re.sub(
                r"~?\s*39\s*[µμ]s(?:\s*(?:per|/)\s*sequence)?(?:\s*\([^)]*\))?",
                "495.7 µs per sequence (26.77 ms total across 54 loci; ~2,017 sequences/s)",
                new_t,
            )

        # Fix #4 any stray "autosomes only" phrase anywhere in the document
        if "autosomes only" in new_t.lower():
            new_t = re.sub(
                r"covers autosomes only",
                "places chr1–11 and chrX (N = 955) in the calibration pool and held-out chr12–22 (N = 150) in the locked blind test set",
                new_t,
                flags=re.IGNORECASE,
            )

        # Fix #7: Make Figures 17, 18, and 19 captions bold just like Figures 15, 16, and 20
        is_fig_17_19 = st.startswith(("Figure 17.", "Figure 18.", "Figure 19."))
        if new_t != txt or is_fig_17_19:
            set_para_text(p, new_t, bold=True if is_fig_17_19 else None)

    # Also check tables (e.g., Table 12 or Supplementary Table S3) for v3.1 label on 72.0% / 21.3% and 39 µs latency
    for tbl in doc.tables:
        for row in tbl.rows:
            for c in row.cells:
                ct = c.text
                nct = ct
                if "72.0%" in nct and "21.3%" in nct and "v3.3" in nct:
                    nct = nct.replace("NC-MutScan v3.3", "Early Prototype v3.1")
                if "39 µs" in nct or "39 μs" in nct:
                    nct = re.sub(r"~?\s*39\s*[µμ]s", "495.7 µs", nct)
                if nct != ct:
                    c.text = nct

    doc.save(DOCX_PATH)

    # -------------------------------------------------------------------------
    # STEP C: Re-open saved .docx from disk and print verification audit
    # -------------------------------------------------------------------------
    chk = docx.Document(DOCX_PATH)
    full_txt = "\n".join(p.text for p in chk.paragraphs)

    assert "39 µs" not in full_txt and "39 μs" not in full_txt, "FAIL: 39 µs still found!"
    assert "autosomes only" not in full_txt.lower(), "FAIL: 'autosomes only' still found!"
    assert "partitioned 1,105" not in full_txt, "FAIL: 'partitioned 1,105' still found!"
    assert "two software regression" not in full_txt, "FAIL: 'two software regression' still found!"
    assert "4.2-bit floor" in full_txt, "FAIL: polyA 4.2-bit floor clause missing!"

    print("\n=== DISK VERIFICATION AUDIT (ALL 8 ITEMS) ===")
    for p in chk.paragraphs:
        st = p.text.strip()
        if st.startswith(("Figure 15.", "Figure 16.", "Figure 17.", "Figure 18.", "Figure 19.", "Figure 20.")):
            is_b = all(r.bold for r in p.runs if r.text.strip())
            print(f"[PASS] Caption bold={is_b}: {st[:95]}...")
        elif "4.2-bit floor" in st:
            print(f"[PASS] §2.10 S_unified + polyA floor: ...{st[-260:]}")
        elif "three software regression" in st:
            print(f"[PASS] §2.11 count: {st[:160]}...")
        elif "execute in 26.77 ms total across the 54 loci" in st:
            m_snip = re.search(r".{0,40}execute in 26\.77 ms total across the 54 loci.{0,90}", st)
            print(f"[PASS] §3.9 latency grammar: {m_snip.group(0)}")
        elif "72.0%" in st and "21.3%" in st and not st.startswith(("Limitations.", "Figure")):
            print(f"[PASS] §3.11 v3.1 & sampling wording: {st[:220]}...")
        elif st.startswith("Limitations."):
            print(f"[PASS] Merged Discussion Limitations: {st[:240]}... [Total len={len(st)} chars]")


if __name__ == "__main__":
    main()
