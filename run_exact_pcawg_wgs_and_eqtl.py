"""
run_exact_pcawg_wgs_and_eqtl.py
Inspects exact headers of pcawg_specimen_histology_August2016_v9 and _donor,
then computes:
  1. PCAWG Skin-Melanoma (MELA-AU) WGS donor counts & UV dipyrimidine signature.
  2. Histology-adjusted & within-B-cell-lymphoma PCAWG eQTL OLS (BCL6, BCL2, PAX5, WDR74).
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


def build_pcawg_maps():
    sample_to_donor = {}
    donor_to_hist = {}
    for fn in ("pcawg_specimen_histology_August2016_v9", "pcawg_specimen_histology_August2016_v9_donor"):
        if not os.path.exists(fn):
            continue
        with open_text_auto(fn) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            print(f"  [{fn}] headers: {hdr}")
            # Find donor column and histology column dynamically
            d_idx = next((i for i, h in enumerate(hdr) if "donor" in h.lower()), 0)
            h_idx = next((i for i, h in enumerate(hdr) if "histology_abbreviation" in h.lower() or "histology" in h.lower() or "project" in h.lower()), None)
            for line in f:
                pts = line.rstrip("\n").split("\t")
                if len(pts) <= d_idx:
                    continue
                did = pts[d_idx].strip()
                if not did:
                    continue
                if h_idx is not None and len(pts) > h_idx and pts[h_idx].strip():
                    donor_to_hist[did] = pts[h_idx].strip()
                for val in pts:
                    v = val.strip()
                    if v and (v.startswith("DO") or v.startswith("SP") or v.startswith("SA") or "-" in v):
                        sample_to_donor[v] = did
    return sample_to_donor, donor_to_hist


def main():
    print("=== 1. PCAWG Histology Labels & Melanoma WGS Check ===")
    sample_to_donor, donor_to_hist = build_pcawg_maps()
    hist_counts = {}
    for h in donor_to_hist.values():
        hist_counts[h] = hist_counts.get(h, 0) + 1
    mela_keys = [h for h in hist_counts if "mela" in h.lower() or "skin" in h.lower()]
    print("  Melanoma histology keys:", {k: hist_counts[k] for k in mela_keys})

    loci = {
        "TERT (C228T/C250T)": ("5", 1295150, 1295300),
        "AC040160.1":         ("16", 67260850, 67261150),
        "AL358472.4":         ("1", 153963050, 153963400),
        "TCF3":               ("19", 1652850, 1653200),
    }
    hits = {k: [] for k in loci}
    mela_donors_seen = set()
    mela_snvs = 0

    pcawg_targets = [
        ("BCL6",  "ENSG00000113916", "3",  187463000, 1500),
        ("BCL2",  "ENSG00000171791", "18", 60985500,  1500),
        ("PAX5",  "ENSG00000196092", "9",  37034000,  1500),
        ("WDR74", "ENSG00000149541", "11", 62609500,  1500),
    ]
    tmb_map = {}
    mut_donors = {t[0]: set() for t in pcawg_targets}

    with open_text_auto("October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS") as f:
        _ = f.readline()
        for line in f:
            pts = line.rstrip("\n").split("\t")
            sid = pts[0].strip()
            did = sample_to_donor.get(sid, sid)
            tmb_map[did] = tmb_map.get(did, 0) + 1
            ch = pts[1].replace("chr", "")
            pos = int(pts[2])
            hist = donor_to_hist.get(did, "")

            if hist in mela_keys:
                mela_donors_seen.add(did)
                mela_snvs += 1
                if ch in ("1", "5", "16", "19"):
                    mut_str = f"{pts[4]}>{pts[5]}"
                    for name, (lc, l0, l1) in loci.items():
                        if ch == lc and l0 <= pos <= l1:
                            hits[name].append((did, pos, mut_str))

            for gene, _, tc, tpos, tol in pcawg_targets:
                if ch == tc and abs(pos - tpos) <= tol:
                    mut_donors[gene].add(did)

    print(f"  PCAWG Melanoma WGS donors in MAF: {len(mela_donors_seen)} ({mela_snvs:,} SNVs)")
    for name, evs in hits.items():
        ud = sorted({e[0] for e in evs})
        ct_ga = sum(1 for e in evs if e[2] in ("C>T", "G>A"))
        positions = sorted({e[1] for e in evs})
        print(f"    {name:20s}: {len(ud):2d}/{len(mela_donors_seen)} donors ({len(ud)/max(1,len(mela_donors_seen))*100:5.1f}%) | "
              f"SNVs={len(evs):2d} | C>T/G>A={ct_ga}/{len(evs)} | pos={positions[:4]}")

    ext_dir = os.path.join("results", "external_pcawg")
    ensg_lookup = {t[1]: t[0] for t in pcawg_targets}
    sym_lookup = {t[0]: t[0] for t in pcawg_targets}

    expr_map = {t[0]: {} for t in pcawg_targets}
    with open_text_auto(os.path.join(ext_dir, "pcawg_rnaseq_fpkm_uq.tsv.gz")) as f:
        cols = [sample_to_donor.get(c.strip(), c.strip()) for c in f.readline().rstrip("\n").split("\t")[1:]]
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
        cols = [sample_to_donor.get(c.strip(), c.strip()) for c in f.readline().rstrip("\n").split("\t")[1:]]
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
                d = sample_to_donor.get(pts[0].strip(), pts[0].strip())
                try:
                    pur_map[d] = float(pts[p_i])
                except ValueError:
                    pass

    matched_donors = set(expr_map["BCL6"].keys()) & set(tmb_map.keys())
    print(f"\n=== 2. Histology-Adjusted & Within-Lymphoma PCAWG eQTL OLS (N={len(matched_donors)} matched) ===")
    for mode in ("Histology-Adjusted (All Matched Cohorts)", "Within B-Cell Lymphoma Only (Lymph-BNHL + Lymph-CLL)"):
        print(f"\n  --- {mode} ---")
        pvals_mode = []
        rows_mode = []
        for gene, _, _, _, _ in pcawg_targets:
            mut_hists = {donor_to_hist.get(d, "Unknown") for d in (mut_donors[gene] & matched_donors)}
            if "Lymphoma" in mode:
                cohort_d = [d for d in matched_donors if "lymph" in donor_to_hist.get(d, "").lower()]
            else:
                cohort_d = [d for d in matched_donors if donor_to_hist.get(d, "Unknown") in mut_hists]

            valid_d = [
                d for d in cohort_d
                if d in expr_map[gene] and d in cnv_map[gene] and d in pur_map and d in tmb_map
                and not math.isnan(pur_map[d]) and not math.isnan(cnv_map[gene][d])
            ]
            hists = sorted({donor_to_hist.get(d, "Unknown") for d in valid_d})
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
                ] + [1.0 if donor_to_hist.get(d) == h else 0.0 for h in hists[1:]]
                X_rows.append(row)
                y_rows.append(expr_map[gene][d])

            if len(X_rows) < 5:
                print(f"    {gene:6s}: insufficient rows ({len(X_rows)}), hists={hists}")
                continue
            b, lo, hi, pv, df = ols_fit(np.array(X_rows), np.array(y_rows), col=1)
            pvals_mode.append(pv)
            rows_mode.append((gene, len(valid_d) - n_mut, n_mut, b, lo, hi, pv))

        qvals = bh_fdr(pvals_mode)
        for (gene, n_wt, n_mut, b, lo, hi, pv), qv in zip(rows_mode, qvals):
            print(f"    {gene:6s} | WT={n_wt:3d} MUT={n_mut:2d} | beta={b:+.3f} [{lo:+.3f}, {hi:+.3f}] | p={pv:.3e} | qBH={qv:.3e}")


if __name__ == "__main__":
    main()
