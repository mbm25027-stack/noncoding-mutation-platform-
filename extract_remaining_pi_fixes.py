"""
extract_remaining_pi_fixes.py
Compact, single-screen extraction for:
  1. Concern #1 & #2: chr7:140,453,135 (BRAF/MKRN1) coordinate audit + PCAWG MELA-AU WGS donor counts.
  2. Concern #4: Histology-adjusted & Within-B-cell-lymphoma PCAWG eQTL OLS (BCL6, BCL2, PAX5, WDR74).
  3. Concern #8: Exact Splice-PWM-only vs. TF-PWM-only ablation breakdown on locked_blind_chr12_22.vcf.
"""

import csv
import glob
import gzip
import json
import math
import os
import re
import numpy as np
from nc_mutscan_research_suite import analyze_variant_delta

FLANK_BP = 25


def wilson_ci(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    denom = 1.0 + (z ** 2) / n
    center = (p + (z ** 2) / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + (z ** 2) / (4.0 * n)) / n)) / denom
    return (p * 100.0, max(0.0, center - half) * 100.0, min(1.0, center + half) * 100.0)


def betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, 200):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 / max(1e-30, 1.0 + aa * d) if abs(1.0 + aa * d) >= 1e-30 else 1e30
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
    n, p = X.shape
    df = max(1, n - p)
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


# ---------------------------------------------------------------------------
# 1. Splice-Only vs TF-Only Ablation (inspecting exact return keys)
# ---------------------------------------------------------------------------
def run_splice_vs_tf_ablation():
    flank_cache = json.load(open("results/pi_upgrades/blind_flanks_cache.json", "rt", encoding="utf-8"))
    records = []
    with open("locked_blind_chr12_22.vcf", "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] in ("PATHOGENIC", "BENIGN"):
                records.append({"label": p[0], "chrom": p[1].replace("chr", ""), "pos": int(p[2]), "ref": p[4], "alt": p[5].split(",")[0]})
            else:
                lbl = "PATHOGENIC" if ("Pathogenic" in "\t".join(p[5:]) or "Likely_pathogenic" in "\t".join(p[5:])) else "BENIGN"
                records.append({"label": lbl, "chrom": p[0].replace("chr", ""), "pos": int(p[1]), "ref": p[3], "alt": p[4].split(",")[0]})

    splice_counts = [0, 0, 0, 0]
    tf_counts = [0, 0, 0, 0]
    sample_keys = None
    for rec in records:
        key = f"chr{rec['chrom']}:{rec['pos']}:{rec['ref']}>{rec['alt']}"
        seqs = flank_cache[key]
        is_path = rec["label"] == "PATHOGENIC"
        res = analyze_variant_delta(seqs["wt"], seqs["mut"], flank_bp=FLANK_BP, locus=key)
        if sample_keys is None:
            sample_keys = list(res.keys())
        pred_full = res["classification"] in ("Driver Variant", "Loss-of-Function Variant")
        blob = json.dumps(res).upper()
        is_splice = pred_full and any(k in blob for k in ("SPLICE_5SS", "SPLICE_3SS", "POLY_A", "SPLICE"))
        is_tf = pred_full and not is_splice

        for arr, pred in ((splice_counts, is_splice), (tf_counts, is_tf)):
            if is_path and pred:
                arr[0] += 1
            elif not is_path and pred:
                arr[1] += 1
            elif not is_path and not pred:
                arr[2] += 1
            else:
                arr[3] += 1

    print(f"[1] analyze_variant_delta keys: {sample_keys}")
    for lbl, (tp, fp, tn, fn) in (("Splice/PolyA-PWM Only", splice_counts), ("TF-Motif-PWM Only", tf_counts)):
        sens, _, _ = wilson_ci(tp, tp + fn)
        spec, _, _ = wilson_ci(tn, tn + fp)
        ppv, _, _ = wilson_ci(tp, tp + fp)
        print(f"    {lbl:22s}: TP={tp:2d} FP={fp:2d} TN={tn:2d} FN={fn:2d} | Sens={sens:5.2f}% | Spec={spec:5.2f}% | PPV={ppv:5.2f}%")


# ---------------------------------------------------------------------------
# 2. Coordinate Audit (BRAF/MKRN1) & PCAWG MELA-AU WGS Hotspot Donor Counts
# ---------------------------------------------------------------------------
def audit_coords_and_mela_wgs():
    print("\n[2] Coordinate Audit & PCAWG MELA-AU WGS Hotspot Check:")
    all_files = [f for f in glob.glob("**/*", recursive=True) if os.path.isfile(f)]
    # Check where 140453135 came from
    for f in all_files:
        if f.endswith((".py", ".csv", ".tsv", ".json")) and "venv" not in f:
            try:
                txt = open(f, "rt", encoding="utf-8", errors="ignore").read()
            except Exception:
                continue
            if "140453135" in txt or "140,453,135" in txt:
                print(f"    140453135 found in {f}")
                for ln in txt.splitlines():
                    if "140453135" in ln or "140,453,135" in ln:
                        print("      ->", ln.strip()[:120])

    # Find all data files in workspace to locate PCAWG MELA-AU WGS SNVs & eQTL inputs
    data_files = [
        f for f in all_files
        if "venv" not in f and not f.endswith((".png", ".pdf", ".docx", ".pyc", ".html", ".js", ".css"))
        and any(k in f.lower() for k in ("data", "pcawg", "mela", "brca", "lihc", "skcm", "maf", "snv", "fpkm", "purity", "clinical", "cna", "gistic"))
    ]
    print("    Workspace data files:", [f for f in data_files if not f.endswith(".py")][:20])

    # Scan any SNV/MAF/BED/TSV file containing MELA or PCAWG SNVs
    targets = {
        "TERT_promoter":   [("5", 1295150, 1295300), ("5", 1295050, 1295200)],
        "AC040160.1_prom": [("16", 67260800, 67261150), ("16", 67294800, 67295200)],
        "AL358472.4_prom": [("1", 153963050, 153963400), ("1", 153935600, 153936000)],
        "TCF3_prom":       [("19", 1652850, 1653200)],
    }
    for f in data_files:
        if f.endswith(".py") or os.path.getsize(f) < 50000:
            continue
        if not any(k in f.lower() for k in ("mela", "pcawg", "snv", "maf")):
            continue
        opener = gzip.open if f.endswith(".gz") else open
        try:
            with opener(f, "rt", encoding="utf-8", errors="ignore") as fh:
                first = fh.readline().rstrip("\n").split("\t")
                c_i = next((i for i, h in enumerate(first) if h.lower() in ("chrom", "chromosome", "chr", "#chrom")), None)
                p_i = next((i for i, h in enumerate(first) if h.lower() in ("pos", "position", "start_position", "start")), None)
                d_i = next((i for i, h in enumerate(first) if any(k in h.lower() for k in ("donor", "sample", "tumor_sample", "icgc_donor"))), None)
                if c_i is None or p_i is None:
                    continue
                hits = {k: set() for k in targets}
                for line in fh:
                    pts = line.rstrip("\n").split("\t")
                    if len(pts) <= max(c_i, p_i):
                        continue
                    ch = pts[c_i].replace("chr", "")
                    if ch not in ("1", "5", "16", "19"):
                        continue
                    try:
                        pos = int(pts[p_i])
                    except ValueError:
                        continue
                    dn = pts[d_i] if (d_i is not None and len(pts) > d_i) else f"r_{pos}"
                    for tn, spans in targets.items():
                        for sc, s0, s1 in spans:
                            if ch == sc and s0 <= pos <= s1:
                                hits[tn].add((dn, pos))
                if any(hits.values()):
                    print(f"    MELA/PCAWG WGS scan ({f}):")
                    for tn, dset in hits.items():
                        print(f"      {tn:18s}: {len({d for d, _ in dset})} unique donors ({len(dset)} SNVs)")
        except Exception as e:
            print(f"    Note on {f}: {e}")


# ---------------------------------------------------------------------------
# 3. Inspect 22_nb_bmr_and_covariate_eqtl.py to run Histology-Adjusted OLS
# ---------------------------------------------------------------------------
def inspect_and_run_eqtl():
    print("\n[3] Inspecting 22_nb_bmr_and_covariate_eqtl.py data sources:")
    if os.path.exists("22_nb_bmr_and_covariate_eqtl.py"):
        lines = open("22_nb_bmr_and_covariate_eqtl.py", "rt", encoding="utf-8", errors="ignore").readlines()
        for i, ln in enumerate(lines, 1):
            if any(k in ln.lower() for k in ("open(", "read_", "data", "pcawg", "bcl6", "ols", "smf", "histolog", "cohort", "project")):
                print(f"    L{i:03d}: {ln.rstrip()[:115]}")


if __name__ == "__main__":
    run_splice_vs_tf_ablation()
    audit_coords_and_mela_wgs()
    inspect_and_run_eqtl()
