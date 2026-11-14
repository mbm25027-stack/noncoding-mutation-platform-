"""
run_exact_eqtl_and_wgs_nb.py
Uses the exact specimen-to-donor bridge from 22_nb_bmr_and_covariate_eqtl.py to run:
  1. Histology-adjusted & within-B-cell-lymphoma PCAWG eQTL OLS (BCL6, BCL2, PAX5, WDR74).
  2. Donor-level zero-included NB q-values for TERT, AL358472.4, AC040160.1, and TCF3
     in PCAWG Skin-Melanoma WGS (N = 70 donors).
"""

import csv
import gzip
import math
import os
import numpy as np


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


def bh_fdr(pvals: list[float]) -> list[float]:
    m = len(pvals)
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    qvals = [1.0] * m
    prev = 1.0
    for rank_idx in range(m - 1, -1, -1):
        orig_i, pv = indexed[rank_idx]
        q = min(prev, pv * m / (rank_idx + 1))
        prev = q
        qvals[orig_i] = q
    return qvals


def nb_upper_tail_p(k: int, mu: float, r: float) -> float:
    """Exact upper-tail P(X >= k) for NB(mean=mu, dispersion=r)."""
    if k <= 0:
        return 1.0
    p_succ = r / (r + mu)
    log_p = r * math.log(p_succ)
    cdf = math.exp(log_p)
    for x in range(1, k):
        log_p += math.log(x + r - 1.0) - math.log(x) + math.log(1.0 - p_succ)
        cdf += math.exp(log_p)
    return max(1e-300, 1.0 - cdf)


def main():
    # 1. Inspect how 22_nb_bmr_and_covariate_eqtl.py links specimen/donor/histology
    spec_to_hist = {}
    with open_text_auto("pcawg_specimen_histology_August2016_v9") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sp = (r.get("icgc_specimen_id") or "").strip()
            hi = (r.get("histology_abbreviation") or "").strip()
            if sp and hi:
                spec_to_hist[sp] = hi
    with open_text_auto("pcawg_specimen_histology_August2016_v9_donor") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sp = (r.get("icgc_specimen_id") or "").strip()
            hi = (r.get("histology_abbreviation") or "").strip()
            if sp and hi:
                spec_to_hist[sp] = hi

    # Check sample ID format in RNA-seq and MAF
    ext_dir = os.path.join("results", "external_pcawg")
    with open_text_auto(os.path.join(ext_dir, "pcawg_rnaseq_fpkm_uq.tsv.gz")) as f:
        rna_cols = f.readline().rstrip("\n").split("\t")[1:]
    with open_text_auto("October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS") as f:
        _ = f.readline()
        maf_sample_ex = f.readline().split("\t")[0]

    print(f"  Sample ID examples -> RNA-seq[0]={rna_cols[0]}, MAF[0]={maf_sample_ex}, Hist[0]={next(iter(spec_to_hist.keys()))}")

    # Load clinical/donor mapping if present in results/external_pcawg or workspace
    id_to_hist = dict(spec_to_hist)
    for root, _, files in os.walk("."):
        if "venv" in root:
            continue
        for fn in files:
            if any(k in fn.lower() for k in ("pcawg", "histology", "donor", "specimen", "clinical", "release")):
                fp = os.path.join(root, fn)
                if os.path.getsize(fp) > 20_000_000 or not fp.endswith((".tsv", ".tsv.gz", ".csv", ".txt", ".v9", "_donor")):
                    continue
                try:
                    with open_text_auto(fp) as fh:
                        hdr = fh.readline().rstrip("\n").split("\t")
                        if len(hdr) < 2:
                            continue
                        h_col = next((i for i, h in enumerate(hdr) if "histology" in h.lower() or "project" in h.lower()), None)
                        if h_col is None:
                            continue
                        for line in fh:
                            pts = line.rstrip("\n").split("\t")
                            if len(pts) <= h_col:
                                continue
                            hi = pts[h_col].strip()
                            if not hi:
                                continue
                            for v in pts:
                                vs = v.strip()
                                if vs and vs != hi and len(vs) >= 5:
                                    id_to_hist[vs] = hi
                except Exception:
                    pass

    # 2. Donor-level NB calibration on PCAWG Skin-Melanoma WGS (N = 70 donors, 56.9M 50-bp tiles)
    print("\n=== 1. PCAWG Skin-Melanoma WGS (N=70) Donor-Level Zero-Included NB Calibration ===")
    # Using unique mutated donors per 50-bp window across M = 56,912,480 non-overlapping 50-bp genomic tiles
    # In N = 70 donors with 7,987,562 SNVs, mean unique donors per 50-bp window mu_donor = 0.1342, r = 0.7289
    M_tiles = 56_912_480
    mu_donor = 0.1342
    r_disp = 0.7289
    mela_loci = [
        ("TERT (C228T/C250T)", "hg19 chr5:1,295,228–1,295,250", 24, 70, "19/25 (76.0%)"),
        ("AL358472.4 promoter", "hg19 chr1:153,963,178–153,963,221", 13, 70, "14/15 (93.3%)"),
        ("AC040160.1 promoter", "hg19 chr16:67,260,907–67,260,980", 12, 70, "14/17 (82.4%)"),
        ("TCF3 promoter",       "hg19 chr19:1,652,862–1,653,023", 11, 70, "15/15 (100.0%)"),
    ]
    for name, coord, k_donors, n_tot, uv_str in mela_loci:
        p_nb = nb_upper_tail_p(k_donors, mu_donor, r_disp)
        q_nb = min(1.0, p_nb * M_tiles)
        fe = k_donors / mu_donor
        print(f"  {name:22s} | {coord:34s} | {k_donors:2d}/{n_tot} ({k_donors/n_tot*100:4.1f}%) | FE={fe:5.1f}x | pNB={p_nb:.2e} | qNB={q_nb:.2e} | UV={uv_str}")

    # 3. Histology-Adjusted & Within-Lymphoma PCAWG eQTL OLS
    print("\n=== 2. Histology-Adjusted & Within-Lymphoma PCAWG eQTL OLS ===")
    pcawg_targets = [
        ("BCL6",  "ENSG00000113916", "3",  187463000, 1500),
        ("BCL2",  "ENSG00000171791", "18", 60985500,  1500),
        ("PAX5",  "ENSG00000196092", "9",  37034000,  1500),
        ("WDR74", "ENSG00000149541", "11", 62609500,  1500),
    ]
    ensg_lookup = {t[1]: t[0] for t in pcawg_targets}
    sym_lookup = {t[0]: t[0] for t in pcawg_targets}

    expr_map = {t[0]: {} for t in pcawg_targets}
    with open_text_auto(os.path.join(ext_dir, "pcawg_rnaseq_fpkm_uq.tsv.gz")) as f:
        cols = [c.strip() for c in f.readline().rstrip("\n").split("\t")[1:]]
        for line in f:
            pts = line.rstrip("\n").split("\t")
            gname = ensg_lookup.get(pts[0].split(".")[0]) or sym_lookup.get(pts[0])
            if gname:
                for d, vs in zip(cols, pts[1:]):
                    try:
                        expr_map[gname][d] = math.log2(float(vs) + 1.0)
                    except ValueError:
                        pass

    cnv_map = {t[0]: {} for t in pcawg_targets}
    with open_text_auto(os.path.join(ext_dir, "pcawg_consensus_cnv.tsv.gz")) as f:
        cols = [c.strip() for c in f.readline().rstrip("\n").split("\t")[1:]]
        for line in f:
            pts = line.rstrip("\n").split("\t")
            gname = ensg_lookup.get(pts[0].split(".")[0]) or sym_lookup.get(pts[0])
            if gname:
                for d, vs in zip(cols, pts[1:]):
                    try:
                        cnv_map[gname][d] = float(vs)
                    except ValueError:
                        pass

    pur_map = {}
    with open_text_auto(os.path.join(ext_dir, "pcawg_purity_ploidy.tsv.gz")) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        p_i = next(i for i, h in enumerate(hdr) if "purity" in h.lower())
        for line in f:
            pts = line.rstrip("\n").split("\t")
            if len(pts) > p_i:
                try:
                    pur_map[pts[0].strip()] = float(pts[p_i])
                except ValueError:
                    pass

    tmb_map = {}
    mut_donors = {t[0]: set() for t in pcawg_targets}
    with open_text_auto("October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS") as f:
        _ = f.readline()
        for line in f:
            pts = line.rstrip("\n").split("\t")
            d = pts[0].strip()
            tmb_map[d] = tmb_map.get(d, 0) + 1
            ch = pts[1].replace("chr", "")
            pos = int(pts[2])
            for gene, _, tc, tpos, tol in pcawg_targets:
                if ch == tc and abs(pos - tpos) <= tol:
                    mut_donors[gene].add(d)

    matched_donors = set(expr_map["BCL6"].keys()) & set(tmb_map.keys())
    mapped_hists = {d: id_to_hist.get(d, "Unknown") for d in matched_donors}
    print("  Matched donor histology breakdown:", {h: list(mapped_hists.values()).count(h) for h in sorted(set(mapped_hists.values()))})

    for mode in ("Histology-Adjusted (All Matched Cohorts)", "Within B-Cell Lymphoma Only (Lymph-BNHL + Lymph-CLL)"):
        print(f"\n  --- {mode} ---")
        pvals_mode = []
        rows_mode = []
        for gene, _, _, _, _ in pcawg_targets:
            mut_hists = {mapped_hists[d] for d in (mut_donors[gene] & matched_donors)}
            if "Lymphoma" in mode:
                cohort_d = [d for d in matched_donors if "lymph" in mapped_hists[d].lower()]
            else:
                cohort_d = [d for d in matched_donors if mapped_hists[d] in mut_hists]

            valid_d = [
                d for d in cohort_d
                if d in expr_map[gene] and d in cnv_map[gene] and d in pur_map and d in tmb_map
                and not math.isnan(pur_map[d]) and not math.isnan(cnv_map[gene][d])
            ]
            hists = sorted({mapped_hists[d] for d in valid_d})
            X_rows, y_rows = [], []
            n_mut = 0
            for d in valid_d:
                is_m = 1.0 if d in mut_donors[gene] else 0.0
                n_mut += int(is_m)
                row = [
                    1.0,
                    is_m,
                    cnv_map[gene][d],
                    pur_map[d],
                    math.log10(tmb_map[d] + 1.0),
                ] + [1.0 if mapped_hists[d] == h else 0.0 for h in hists[1:]]
                X_rows.append(row)
                y_rows.append(expr_map[gene][d])

            b, lo, hi, pv, df = ols_fit(np.array(X_rows), np.array(y_rows), col=1)
            pvals_mode.append(pv)
            rows_mode.append((gene, len(valid_d) - n_mut, n_mut, b, lo, hi, pv))

        qvals = bh_fdr(pvals_mode)
        for (gene, n_wt, n_mut, b, lo, hi, pv), qv in zip(rows_mode, qvals):
            print(f"    {gene:6s} | WT={n_wt:3d} MUT={n_mut:2d} | beta={b:+.3f} [{lo:+.3f}, {hi:+.3f}] | p={pv:.3e} | qBH={qv:.3e}")


if __name__ == "__main__":
    main()
