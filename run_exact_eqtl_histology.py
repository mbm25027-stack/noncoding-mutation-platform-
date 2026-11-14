"""
run_exact_eqtl_histology.py
Maps SP... <-> DO... across all metadata files, overlaps segment-level CNV calls,
and computes Histology-Adjusted, Within-B-Cell-Lymphoma, and Within-Lymph-BNHL eQTL OLS.
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


def main():
    # 1. Build complete SP/SA/UUID -> DO and DO -> histology_abbreviation maps
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

    ext_dir = os.path.join("results", "external_pcawg")
    # (Gene, ENSG, chrom, prom_center_hg19, prom_tol, gene_start_hg19, gene_end_hg19)
    pcawg_targets = [
        ("BCL6",  "ENSG00000113916", "3",  187463000, 1500, 187439165, 187463515),
        ("BCL2",  "ENSG00000171791", "18", 60985500,  1500, 60790579,  60987361),
        ("PAX5",  "ENSG00000196092", "9",  37034000,  1500, 36833272,  37034476),
        ("WDR74", "ENSG00000149541", "11", 62609500,  1500, 62608568,  62611395),
    ]
    ensg_lookup = {t[1]: t[0] for t in pcawg_targets}
    sym_lookup = {t[0]: t[0] for t in pcawg_targets}

    # 2. Load RNA-seq expression (log2(FPKM-UQ + 1))
    expr_map = {t[0]: {} for t in pcawg_targets}
    with open_text_auto(os.path.join(ext_dir, "pcawg_rnaseq_fpkm_uq.tsv.gz")) as f:
        hdr = [c.strip().strip('"') for c in f.readline().rstrip("\n").split("\t")]
        cols = [id_to_donor.get(c, c) for c in hdr[1:]]
        for line in f:
            pts = line.rstrip("\n").split("\t")
            base = pts[0].strip().strip('"').split(".")[0]
            gname = ensg_lookup.get(base) or sym_lookup.get(base)
            if gname:
                for d, vs in zip(cols, pts[1:]):
                    try:
                        expr_map[gname][d] = math.log2(float(vs.strip('"')) + 1.0)
                    except ValueError:
                        pass

    # 3. Load Segment-Level Consensus CNV overlapping each gene locus
    cnv_map = {t[0]: {} for t in pcawg_targets}
    with open_text_auto(os.path.join(ext_dir, "pcawg_consensus_cnv.tsv.gz")) as f:
        _ = f.readline()  # ['sampleID', 'chr', 'start', 'end', 'value']
        for line in f:
            pts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if len(pts) < 5:
                continue
            ch = pts[1].replace("chr", "")
            if ch not in ("3", "9", "11", "18"):
                continue
            try:
                s0, s1, val = int(pts[2]), int(pts[3]), float(pts[4])
            except ValueError:
                continue
            d = id_to_donor.get(pts[0], pts[0])
            for gene, _, gc, _, _, gs, ge in pcawg_targets:
                if ch == gc and s0 <= ge and s1 >= gs:
                    cnv_map[gene][d] = val

    # 4. Load Purity
    pur_map = {}
    with open_text_auto(os.path.join(ext_dir, "pcawg_purity_ploidy.tsv.gz")) as f:
        hdr = [c.strip().strip('"') for c in f.readline().rstrip("\n").split("\t")]
        p_i = next(i for i, h in enumerate(hdr) if "purity" in h.lower())
        for line in f:
            pts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if len(pts) <= p_i:
                continue
            try:
                val = float(pts[p_i])
            except ValueError:
                continue
            d = id_to_donor.get(pts[0], pts[0])
            pur_map[d] = val

    # 5. Scan WGS SNVs for TMB & promoter mutation status
    tmb_map = {}
    mut_donors = {t[0]: set() for t in pcawg_targets}
    with open_text_auto("October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS") as f:
        _ = f.readline()
        for line in f:
            pts = line.rstrip("\n").split("\t")
            d = id_to_donor.get(pts[0].strip(), pts[0].strip())
            tmb_map[d] = tmb_map.get(d, 0) + 1
            ch = pts[1].replace("chr", "")
            pos = int(pts[2])
            for gene, _, tc, tpos, tol, _, _ in pcawg_targets:
                if ch == tc and abs(pos - tpos) <= tol:
                    mut_donors[gene].add(d)

    matched_donors = set(expr_map["BCL6"].keys()) & set(tmb_map.keys())
    print(f"Matched={len(matched_donors)} | Expr={len(expr_map['BCL6'])} | CNV_BCL6={len(cnv_map['BCL6'])} | Pur={len(pur_map)}")

    os.makedirs("results/clean_submission", exist_ok=True)
    out_rows = []
    modes = [
        ("1. Histology-Adjusted (+ histology_abbreviation fixed effects)", "all_matched"),
        ("2. Within B-Cell Lymphoma Only (Lymph-BNHL + Lymph-CLL)", "lymphoma_both"),
        ("3. Within Lymph-BNHL Alone (Germinal-Center B-NHL Only)", "bnhl_only"),
    ]
    for title, mkey in modes:
        print(f"\n--- {title} ---")
        pvals_mode = []
        rows_mode = []
        for gene, _, _, _, _, _, _ in pcawg_targets:
            mut_hists = {donor_to_hist.get(d, "Unknown") for d in (mut_donors[gene] & matched_donors)}
            if mkey == "lymphoma_both":
                cohort_d = [d for d in matched_donors if donor_to_hist.get(d) in ("Lymph-BNHL", "Lymph-CLL")]
            elif mkey == "bnhl_only":
                cohort_d = [d for d in matched_donors if donor_to_hist.get(d) == "Lymph-BNHL"]
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
                ] + [1.0 if donor_to_hist.get(d) == hi else 0.0 for hi in hists[1:]]
                X_rows.append(row)
                y_rows.append(expr_map[gene][d])

            if len(X_rows) < 5 or n_mut == 0:
                print(f"  {gene:6s}: skipped (N={len(X_rows)}, MUT={n_mut})")
                continue
            b, lo, hi, pv, df = ols_fit(np.array(X_rows), np.array(y_rows), col=1)
            pvals_mode.append(pv)
            rows_mode.append((gene, len(valid_d) - n_mut, n_mut, b, lo, hi, pv, hists))

        qvals = bh_fdr(pvals_mode)
        for (gene, n_wt, n_mut, b, lo, hi, pv, hists), qv in zip(rows_mode, qvals):
            print(f"  {gene:6s} | WT={n_wt:3d} MUT={n_mut:2d} | beta={b:+.3f} [{lo:+.3f}, {hi:+.3f}] | p={pv:.3e} | qBH={qv:.3e} | cohorts={hists}")
            out_rows.append({
                "Model": title, "Gene": gene, "WT_Donors": n_wt, "Mut_Donors": n_mut,
                "Beta_Mut": f"{b:.4f}", "CI_Low": f"{lo:.4f}", "CI_High": f"{hi:.4f}",
                "P_Value": f"{pv:.4e}", "Q_BH": f"{qv:.4e}", "Cohorts": ";".join(hists)
            })

    if out_rows:
        out_csv = "results/clean_submission/pcawg_histology_adjusted_eqtl.csv"
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
            w.writeheader()
            w.writerows(out_rows)
        print(f"\nSaved {out_csv}")


if __name__ == "__main__":
    main()
