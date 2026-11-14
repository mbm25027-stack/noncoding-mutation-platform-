"""
audit_and_compute_pi_v2.py
Completes:
  1. PCAWG BCL2 cis-eQTL adjusted for t(14;18) / BCL2 SVs + WT=134 vs WT=120 reconciliation.
  2. Saturation-mutagenesis benchmark of NC-MutScan v3.3 TF-motif Delta-PWM across
     the TERT core promoter (Kircher et al. 2019 regulatory region).
  3. Full structural index of NC-DriverEngine_v3_1_Day6_Final_Manuscript_SUBMISSION_READY.docx.
"""

import gzip
import json
import math
import os
import re
import urllib.request
import numpy as np
import docx
from nc_mutscan_research_suite import analyze_variant_delta


def open_text_auto(path: str):
    with open(path, "rb") as raw:
        magic = raw.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) >= 1e-30 else 1e-30)
    h = d
    for m in range(1, 200):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 / (1.0 + aa * d if abs(1.0 + aa * d) >= 1e-30 else 1e-30)
        c = 1.0 + aa / (c if abs(c) >= 1e-30 else 1e-30)
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 / (1.0 + aa * d if abs(1.0 + aa * d) >= 1e-30 else 1e-30)
        c = 1.0 + aa / (c if abs(c) >= 1e-30 else 1e-30)
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    return h


def betainc_reg(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * betacf(a, b, x) / a
    return 1.0 - front * betacf(b, a, 1.0 - x) / b


def ols_fit(X: np.ndarray, y: np.ndarray, col: int = 1):
    n = X.shape[0]
    rank = int(np.linalg.matrix_rank(X))
    df = max(1, n - rank)
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    s2 = float(np.sum(resid ** 2)) / df
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = math.sqrt(max(1e-18, float(cov[col, col])))
    b_val = float(beta[col])
    t_stat = abs(b_val) / se
    p_val = betainc_reg(df / 2.0, 0.5, df / (df + t_stat ** 2))
    t975 = 1.96 + (2.37 / df)
    return b_val, b_val - t975 * se, b_val + t975 * se, p_val, df


def check_bcl2_with_t14_18():
    print("=== 1. PCAWG BCL2 cis-eQTL Adjusted for t(14;18) / BCL2 SVs ===")
    id_to_donor = {}
    donor_to_hist = {}
    for root, _, files in os.walk("."):
        if "venv" in root:
            continue
        for fn in files:
            fp = os.path.join(root, fn)
            if os.path.getsize(fp) > 12_000_000:
                continue
            if not any(k in fn.lower() for k in ("pcawg", "donor", "specimen", "clinical", "histology", "whitelist")):
                continue
            try:
                with open_text_auto(fp) as fh:
                    hdr = [h.strip().strip('"') for h in fh.readline().rstrip("\n").split("\t")]
                    h_idx = next((i for i, h in enumerate(hdr) if h == "histology_abbreviation"), None)
                    for line in fh:
                        pts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
                        do_id = next((v for v in pts if v.startswith("DO")), "")
                        if do_id:
                            if h_idx is not None and len(pts) > h_idx and pts[h_idx]:
                                donor_to_hist[do_id] = pts[h_idx]
                            for v in pts:
                                if v:
                                    id_to_donor[v] = do_id
            except Exception:
                pass

    sv_bcl2_donors = set()
    t14_18_donors = set()
    sv_file = "pcawg_consensus_1.6.161116.somatic_svs.xena.donor"
    if os.path.exists(sv_file):
        with open_text_auto(sv_file) as f:
            _ = f.readline()
            for line in f:
                pts = line.rstrip("\n").split("\t")
                if len(pts) < 8:
                    continue
                d = id_to_donor.get(pts[0].strip(), pts[0].strip())
                ch = pts[1].replace("chr", "")
                try:
                    pos = int(pts[2])
                except ValueError:
                    pos = -1
                alt_bnd = pts[5]
                g1, g2 = pts[6].upper(), pts[7].upper()
                has_bcl2 = ("BCL2" in g1 or "BCL2" in g2 or (ch == "18" and 60_500_000 <= pos <= 61_500_000) or "18:60" in alt_bnd or "18:61" in alt_bnd)
                if has_bcl2:
                    sv_bcl2_donors.add(d)
                    if "IGH" in g1 or "IGH" in g2 or "14:" in alt_bnd or ch == "14":
                        t14_18_donors.add(d)

    ext_dir = os.path.join("results", "external_pcawg")
    expr_bcl2 = {}
    with open_text_auto(os.path.join(ext_dir, "pcawg_rnaseq_fpkm_uq.tsv.gz")) as f:
        hdr = [c.strip().strip('"') for c in f.readline().rstrip("\n").split("\t")]
        cols = [id_to_donor.get(c, c) for c in hdr[1:]]
        for line in f:
            pts = line.rstrip("\n").split("\t")
            if pts[0].strip().strip('"').split(".")[0] == "ENSG00000171791":
                for d, vs in zip(cols, pts[1:]):
                    try:
                        expr_bcl2[d] = math.log2(float(vs.strip('"')) + 1.0)
                    except ValueError:
                        pass

    cnv_bcl2 = {}
    with open_text_auto(os.path.join(ext_dir, "pcawg_consensus_cnv.tsv.gz")) as f:
        _ = f.readline()
        for line in f:
            pts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if len(pts) >= 5 and pts[1].replace("chr", "") == "18":
                try:
                    s0, s1, val = int(pts[2]), int(pts[3]), float(pts[4])
                except ValueError:
                    continue
                if s0 <= 60987361 and s1 >= 60790579:
                    cnv_bcl2[id_to_donor.get(pts[0], pts[0])] = val

    pur_map = {}
    with open_text_auto(os.path.join(ext_dir, "pcawg_purity_ploidy.tsv.gz")) as f:
        hdr = [c.strip().strip('"') for c in f.readline().rstrip("\n").split("\t")]
        p_i = next(i for i, h in enumerate(hdr) if "purity" in h.lower())
        for line in f:
            pts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if len(pts) > p_i:
                try:
                    pur_map[id_to_donor.get(pts[0], pts[0])] = float(pts[p_i])
                except ValueError:
                    pass

    tmb_map = {}
    bcl2_snv_donors = set()
    with open_text_auto("October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS") as f:
        _ = f.readline()
        for line in f:
            pts = line.rstrip("\n").split("\t")
            d = id_to_donor.get(pts[0].strip(), pts[0].strip())
            tmb_map[d] = tmb_map.get(d, 0) + 1
            if pts[1].replace("chr", "") == "18" and abs(int(pts[2]) - 60985500) <= 1500:
                bcl2_snv_donors.add(d)

    lymph_rna_only = [d for d in expr_bcl2 if donor_to_hist.get(d) in ("Lymph-BNHL", "Lymph-CLL")]
    lymph_all = [d for d in (set(expr_bcl2) & set(tmb_map)) if donor_to_hist.get(d) in ("Lymph-BNHL", "Lymph-CLL")]
    lymph_cc = [d for d in lymph_all if d in cnv_bcl2 and d in pur_map and not math.isnan(pur_map[d]) and not math.isnan(cnv_bcl2[d])]
    print(f"  Lymphoma RNA-seq donors (including 14 without WGS MAF): total={len(lymph_rna_only)} (WT={sum(d not in bcl2_snv_donors for d in lymph_rna_only)}, MUT={sum(d in bcl2_snv_donors for d in lymph_rna_only)})")
    print(f"  Lymphoma matched WGS+RNA: total={len(lymph_all)} (WT={sum(d not in bcl2_snv_donors for d in lymph_all)}, MUT={sum(d in bcl2_snv_donors for d in lymph_all)})")
    print(f"  Complete-case (CNV+Purity): total={len(lymph_cc)} (WT={sum(d not in bcl2_snv_donors for d in lymph_cc)}, MUT={sum(d in bcl2_snv_donors for d in lymph_cc)})")
    print(f"  BCL2 SV donors in complete-case: any_BCL2_SV={sum(d in sv_bcl2_donors for d in lymph_cc)}, t(14;18)={sum(d in t14_18_donors for d in lymph_cc)}")

    # Model A: With BCL2 SV / t(14;18) covariate
    X_sv, y_sv = [], []
    for d in lymph_cc:
        is_m = 1.0 if d in bcl2_snv_donors else 0.0
        is_sv = 1.0 if d in sv_bcl2_donors else 0.0
        is_bnhl = 1.0 if donor_to_hist.get(d) == "Lymph-BNHL" else 0.0
        X_sv.append([1.0, is_m, is_sv, cnv_bcl2[d], pur_map[d], math.log10(tmb_map[d] + 1.0), is_bnhl])
        y_sv.append(expr_bcl2[d])
    b_m, lo_m, hi_m, p_m, _ = ols_fit(np.array(X_sv), np.array(y_sv), col=1)
    b_s, lo_s, hi_s, p_s, _ = ols_fit(np.array(X_sv), np.array(y_sv), col=2)
    print(f"  OLS + BCL2_SV covariate: SNV beta={b_m:+.3f} [{lo_m:+.3f}, {hi_m:+.3f}] p={p_m:.3e} | BCL2_SV beta={b_s:+.3f} [{lo_s:+.3f}, {hi_s:+.3f}] p={p_s:.3e}")

    # Model B: Strictly SV-negative Lymphoma donors (excluding t(14;18) / BCL2 SVs)
    no_sv = [d for d in lymph_cc if d not in sv_bcl2_donors]
    X_ns, y_ns = [], []
    n_mut_ns = 0
    for d in no_sv:
        is_m = 1.0 if d in bcl2_snv_donors else 0.0
        n_mut_ns += int(is_m)
        is_bnhl = 1.0 if donor_to_hist.get(d) == "Lymph-BNHL" else 0.0
        X_ns.append([1.0, is_m, cnv_bcl2[d], pur_map[d], math.log10(tmb_map[d] + 1.0), is_bnhl])
        y_ns.append(expr_bcl2[d])
    b_n, lo_n, hi_n, p_n, _ = ols_fit(np.array(X_ns), np.array(y_ns), col=1)
    print(f"  SV-Negative Lymphoma Only (WT={len(no_sv)-n_mut_ns}, MUT={n_mut_ns}): SNV beta={b_n:+.3f} [{lo_n:+.3f}, {hi_n:+.3f}] p={p_n:.3e}")


# ---------------------------------------------------------------------------
# 2. TERT Core-Promoter Saturation Mutagenesis Benchmark (Kircher et al. 2019 locus)
# ---------------------------------------------------------------------------
def benchmark_tert_saturation_mutagenesis():
    print("\n=== 2. TERT Core-Promoter Saturation Mutagenesis (Kircher et al. 2019 Locus) ===")
    url_fwd = "https://rest.ensembl.org/sequence/region/human/5:1294970..1295270:1?content-type=application/json"
    req_f = urllib.request.Request(url_fwd, headers={"User-Agent": "NC-DriverEngine-Audit/3.3"})
    with urllib.request.urlopen(req_f, timeout=15) as resp:
        fwd_seq = json.loads(resp.read().decode("utf-8"))["seq"].upper()

    activating_sites = {1295113: ("A", "T"), 1295135: ("A",), 1295046: ("G",)}
    tp = fp = tn = fn = 0
    for idx in range(25, len(fwd_seq) - 25):
        gpos = 1294970 + idx
        ref = fwd_seq[idx]
        for alt in ("A", "C", "G", "T"):
            if alt == ref:
                continue
            wt = fwd_seq[idx - 25: idx + 26]
            mut = fwd_seq[idx - 25: idx] + alt + fwd_seq[idx + 1: idx + 26]
            res = analyze_variant_delta(wt, mut, flank_bp=25, locus=f"chr5:{gpos}:{ref}>{alt}")
            pred_driver = res["classification"] in ("Driver Variant", "Loss-of-Function Variant")
            is_act = gpos in activating_sites and alt in activating_sites[gpos]
            if is_act and pred_driver:
                tp += 1
            elif not is_act and pred_driver:
                fp += 1
            elif not is_act and not pred_driver:
                tn += 1
            else:
                fn += 1
    tot = tp + fp + tn + fn
    print(f"  TERT Saturation Mutagenesis (N={tot} SNVs across 251-bp core promoter):")
    print(f"  Canonical activating ETS-gain SNVs detected: {tp}/{tp+fn} ({tp/max(1,tp+fn)*100:.1f}%) | Non-activating specificity: {tn}/{tn+fp} ({tn/max(1,tn+fp)*100:.2f}%) | FP={fp}")


# ---------------------------------------------------------------------------
# 3. Structural Audit of Clean Original Manuscript .docx
# ---------------------------------------------------------------------------
def dump_original_manuscript_outline():
    src = "NC-DriverEngine_v3_1_Day6_Final_Manuscript_SUBMISSION_READY.docx"
    doc = docx.Document(src)
    os.makedirs("results/clean_submission", exist_ok=True)
    out_txt = "results/clean_submission/original_manuscript_map.txt"
    with open(out_txt, "w", encoding="utf-8") as out:
        out.write(f"TOTAL PARAGRAPHS: {len(doc.paragraphs)} | TOTAL TABLES: {len(doc.tables)}\n\n")
        for i, p in enumerate(doc.paragraphs):
            t = p.text.strip()
            if not t:
                continue
            st = p.style.name if p.style else ""
            if st.startswith("Heading") or re.match(r"^(\d+\.|\d+\.\d+|Abstract|Introduction|Methods|Results|Discussion|Conclusion|Declarations|References|Table \d+|Figure \d+)", t, re.I):
                out.write(f"P{i:03d} [{st}]: {t[:120]}\n")
        out.write("\n--- TABLES ---\n")
        for ti, tbl in enumerate(doc.tables):
            hdr = [c.text.strip().replace("\n", " ") for c in tbl.rows[0].cells]
            out.write(f"T{ti:02d} (rows={len(tbl.rows)}): {hdr}\n")
    print(f"\n=== 3. Structural Map Saved to {out_txt} ===")
    print(open(out_txt, "rt", encoding="utf-8").read())


if __name__ == "__main__":
    check_bcl2_with_t14_18()
    benchmark_tert_saturation_mutagenesis()
    dump_original_manuscript_outline()
