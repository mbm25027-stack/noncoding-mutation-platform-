import csv
import glob
import gzip
import json
import math
import os
import urllib.request


# ==============================================================================
# 1. Pure-Python OLS & Benjamini-Hochberg Helpers (Zero External DLLs)
# ==============================================================================

def _betacf(a: float, b: float, x: float) -> float:
    max_iter, eps, fpmin = 200, 3.0e-14, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def bh_fdr(p_values: list[float]) -> list[float]:
    m = len(p_values)
    if m == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    q_vals = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        orig_idx, p = indexed[rank - 1]
        q = min(prev, p * m / rank)
        q_vals[orig_idx] = q
        prev = q
    return q_vals


def _invert_matrix(A: list[list[float]]) -> list[list[float]]:
    n = len(A)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(A)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(aug[r][col]))
        aug[col], aug[piv] = aug[piv], aug[col]
        div = aug[col][col]
        if abs(div) < 1e-15:
            raise ValueError("Singular design matrix in OLS regression")
        aug[col] = [v / div for v in aug[col]]
        for r in range(n):
            if r != col:
                factor = aug[r][col]
                aug[r] = [aug[r][c] - factor * aug[col][c] for c in range(2 * n)]
    return [row[n:] for row in aug]


def fit_ols(X: list[list[float]], y: list[float]) -> tuple[list[float], list[float]]:
    n, p = len(y), len(X[0])
    XtX = [[sum(X[i][r] * X[i][c] for i in range(n)) for c in range(p)] for r in range(p)]
    Xty = [sum(X[i][r] * y[i] for i in range(n)) for r in range(p)]
    XtX_inv = _invert_matrix(XtX)
    beta = [sum(XtX_inv[r][c] * Xty[c] for c in range(p)) for r in range(p)]
    rss = sum((y[i] - sum(X[i][c] * beta[c] for c in range(p))) ** 2 for i in range(n))
    df = n - p
    sigma2 = rss / df
    se = [math.sqrt(max(0.0, sigma2 * XtX_inv[c][c])) for c in range(p)]
    t_vals = [beta[c] / se[c] if se[c] > 0 else 0.0 for c in range(p)]
    p_vals = [betainc(0.5 * df, 0.5, df / (df + tv * tv)) for tv in t_vals]
    return beta, p_vals


def open_text_auto(path: str):
    with open(path, "rb") as raw:
        magic = raw.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


# ==============================================================================
# 2. Part A: TCGA-SKCM Covariate-Adjusted eQTL (Using Cached Xena + cBioPortal)
# ==============================================================================

def fetch_cbioportal_skcm_tert_status() -> dict[str, int]:
    """
    Queries open-access cBioPortal REST API (skcm_tcga_pub, TCGA Cell 2015)
    for orthogonal TERT promoter mutation calls (C228T / C250T / CC>TT).
    Returns dict mapping 12-char TCGA patient barcode -> 1 (mutated) or 0 (WT).
    """
    status_by_pat: dict[str, int] = {}
    url = "https://www.cbioportal.org/api/studies/skcm_tcga_pub/clinical-data?clinicalDataType=SAMPLE&pageSize=10000"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=45) as resp:
            records = json.loads(resp.read().decode("utf-8"))
        for rec in records:
            attr = (rec.get("clinicalAttributeId") or "").upper()
            if "TERT" in attr:
                pat = (rec.get("patientId") or "")[:12]
                val = (rec.get("value") or "").strip().lower()
                if pat.startswith("TCGA-"):
                    if val in ("wt", "wildtype", "wild-type", "no", "none", "0", "false"):
                        status_by_pat[pat] = 0
                    elif val not in ("", "na", "nan", "unknown", "not profiled"):
                        status_by_pat[pat] = 1
    except Exception as exc:
        print(f"  [Note] cBioPortal API unavailable ({exc}); relying on local MAF calls.")
    return status_by_pat


def run_tcga_skcm_eqtl() -> list[dict]:
    ext_dir = os.path.join("results", "external_tcga")
    expr_gz = os.path.join(ext_dir, "tcga_skcm_hiseqv2.gz")
    cnv_gz = os.path.join(ext_dir, "tcga_skcm_gistic2.gz")
    clin_txt = os.path.join(ext_dir, "tcga_pancan_purity.txt")
    mut_gz = os.path.join("results", "all_mutations.csv.gz")

    if not (os.path.exists(expr_gz) and os.path.exists(mut_gz)):
        print("[Skip Part A] Cached TCGA-SKCM files not found.")
        return []

    purity_by_pat: dict[str, float] = {}
    if os.path.exists(clin_txt):
        with open_text_auto(clin_txt) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                samp = (row.get("sampleID") or row.get("bcr_patient_barcode") or "").strip()
                if len(samp) < 12:
                    continue
                pat = samp[:12]
                raw_p = row.get("purity") or row.get("percent_tumor_cells") or row.get("percent_tumor_nuclei") or ""
                try:
                    pv = float(raw_p)
                    if pv > 1.0:
                        pv *= 0.01
                    if 0.0 < pv <= 1.0 and pat not in purity_by_pat:
                        purity_by_pat[pat] = pv
                except ValueError:
                    continue

    uv_counts: dict[str, int] = {}
    skcm_pats: set[str] = set()
    wxs_vars: list[tuple[str, str, int]] = []
    with gzip.open(mut_gz, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pat = (row.get("patient") or "")[:12]
            if not pat:
                continue
            skcm_pats.add(pat)
            ref = row.get("Reference_Allele", "")
            alt = row.get("Tumor_Seq_Allele2", "")
            if (ref == "C" and alt == "T") or (ref == "G" and alt == "A"):
                uv_counts[pat] = uv_counts.get(pat, 0) + 1
            chrom = (row.get("Chromosome") or "").replace("chr", "").strip()
            try:
                pos = int(float(row.get("Start_Position", 0)))
                wxs_vars.append((pat, chrom, pos))
            except ValueError:
                continue

    pur_sorted = sorted(purity_by_pat.values())
    med_pur = pur_sorted[len(pur_sorted) // 2] if pur_sorted else 0.75
    print(f"[TCGA-SKCM] Loaded {len(skcm_pats)} patients ({len(purity_by_pat)} with purity; median={med_pur:.2f})")

    targets = [
        ("TERT", "5", 1295113),
        ("RPL13A", "19", 49990694),
        ("DPH3", "3", 16306505),
        ("PES1", "22", 30956798),
        ("MRPS31", "13", 41345136),
        ("ALDOA", "16", 30064190),
    ]
    target_set = {t[0] for t in targets}

    def read_xena_sub(path: str) -> dict[str, dict[str, float]]:
        out: dict[str, dict[str, float]] = {}
        with open_text_auto(path) as f:
            cols = f.readline().rstrip("\n").split("\t")[1:]
            for line in f:
                idx = line.find("\t")
                if idx == -1:
                    continue
                g = line[:idx].split("|")[0].strip()
                if g not in target_set:
                    continue
                vals = line.rstrip("\n").split("\t")[1:]
                pmap: dict[str, float] = {}
                for s, vs in zip(cols, vals):
                    if len(s) >= 15 and s[13:15] in ("01", "06"):
                        p = s[:12]
                        try:
                            fv = float(vs)
                            if not math.isnan(fv) and p not in pmap:
                                pmap[p] = fv
                        except ValueError:
                            continue
                out[g] = pmap
        return out

    expr_map = read_xena_sub(expr_gz)
    cnv_map = read_xena_sub(cnv_gz) if os.path.exists(cnv_gz) else {}
    tert_cbio = fetch_cbioportal_skcm_tert_status()
    if tert_cbio:
        print(f"[TCGA-SKCM] Retrieved orthogonal TERT promoter status for {len(tert_cbio)} patients from cBioPortal")

    results = []
    for gene, chrom, center in targets:
        if gene not in expr_map:
            continue
        if gene == "TERT" and tert_cbio:
            cohort_iter = [p for p in skcm_pats if p in tert_cbio and p in expr_map[gene]]
            mut_set = {p for p in cohort_iter if tert_cbio[p] == 1}
        else:
            cohort_iter = [p for p in skcm_pats if p in expr_map[gene]]
            mut_set = {p for p, c, pos in wxs_vars if c == chrom and abs(pos - center) <= 25}

        X, y, wt_v, mut_v = [], [], [], []
        for p in cohort_iter:
            yv = expr_map[gene][p]
            is_m = 1.0 if p in mut_set else 0.0
            cnv_v = cnv_map.get(gene, {}).get(p, 0.0)
            pur_v = purity_by_pat.get(p, med_pur)
            uv_v = math.log10(uv_counts.get(p, 0) + 1.0)
            X.append([1.0, is_m, cnv_v, pur_v, uv_v])
            y.append(yv)
            (mut_v if is_m == 1.0 else wt_v).append(yv)

        if len(mut_v) < 3:
            continue

        beta, p_vals = fit_ols(X, y)
        n0, n1 = len(wt_v), len(mut_v)
        m0, m1 = sum(wt_v) / n0, sum(mut_v) / n1
        v0 = sum((v - m0) ** 2 for v in wt_v) / (n0 - 1)
        v1 = sum((v - m1) ** 2 for v in mut_v) / (n1 - 1)
        sp2 = ((n0 - 1) * v0 + (n1 - 1) * v1) / (n0 + n1 - 2)
        t_u = (m1 - m0) / math.sqrt(sp2 * (1.0 / n0 + 1.0 / n1))
        df_u = n0 + n1 - 2
        p_u = betainc(0.5 * df_u, 0.5, df_u / (df_u + t_u * t_u))

        results.append({
            "Dataset": "TCGA-SKCM",
            "Gene": gene,
            "Locus_GRCh38": f"chr{chrom}:{center}",
            "WT_Donors": n0,
            "Mut_Donors": n1,
            "Mean_WT_Expr": round(m0, 3),
            "Mean_Mut_Expr": round(m1, 3),
            "Unadj_Log2FC": round(m1 - m0, 3),
            "Unadj_P": p_u,
            "CovAdj_Beta_Mut": round(beta[1], 3),
            "CovAdj_P": p_vals[1],
        })
    return results


# ==============================================================================
# 3. Part B: Open-Access PCAWG WGS + RNA-seq Covariate eQTL (Cached Files)
# ==============================================================================

def run_pcawg_nonus_eqtl() -> list[dict]:
    ext_dir = os.path.join("results", "external_pcawg")
    expr_path = os.path.join(ext_dir, "pcawg_rnaseq_fpkm_uq.tsv.gz")
    cnv_path = os.path.join(ext_dir, "pcawg_consensus_cnv.tsv.gz")
    pur_path = os.path.join(ext_dir, "pcawg_purity_ploidy.tsv.gz")

    maf_hits = [
        p for p in glob.glob("**/*whitelist*.nonUS*", recursive=True)
        if "\\venv\\" not in p and "/venv/" not in p
    ]
    if not (maf_hits and os.path.exists(expr_path)):
        print("[Skip Part B] Local PCAWG nonUS MAF or cached RNA-seq not found.")
        return []

    maf_path = maf_hits[0]

    sp_to_donor: dict[str, str] = {}
    donor_hist: dict[str, str] = {}
    if os.path.exists("pcawg_specimen_histology_August2016_v9") and os.path.exists("pcawg_specimen_histology_August2016_v9_donor"):
        with open_text_auto("pcawg_specimen_histology_August2016_v9") as f_sp, open_text_auto("pcawg_specimen_histology_August2016_v9_donor") as f_do:
            for rsp, rdo in zip(csv.DictReader(f_sp, delimiter="\t"), csv.DictReader(f_do, delimiter="\t")):
                sp = (rsp.get("icgc_specimen_id") or "").strip()
                do = (rdo.get("icgc_specimen_id") or "").strip()
                hi = (rdo.get("histology_abbreviation") or "").strip()
                if sp and do:
                    sp_to_donor[sp] = do
                    sp_to_donor[do] = do
                if do and hi:
                    donor_hist[do] = hi

    donor_purity: dict[str, float] = {}
    if os.path.exists(pur_path):
        with open_text_auto(pur_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                samp = (row.get("samplename") or "").strip()
                do = sp_to_donor.get(samp, samp if samp.startswith("DO") else "")
                try:
                    pv = float(row.get("purity", ""))
                    if do and 0.0 < pv <= 1.0 and do not in donor_purity:
                        donor_purity[do] = pv
                except ValueError:
                    continue

    with open_text_auto(expr_path) as f:
        rna_donors = set(f.readline().rstrip("\n").split("\t")[1:])

    tmb_counts: dict[str, int] = {}
    wgs_vars: list[tuple[str, str, int]] = []
    with open_text_auto(maf_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        sample_col = next(c for c in reader.fieldnames if any(k in c.lower() for k in ("sample", "donor", "patient")))
        chr_col = next(c for c in reader.fieldnames if c.lower() in ("chr", "chrom", "chromosome", "#chrom"))
        pos_col = next(c for c in reader.fieldnames if "start" in c.lower())
        for row in reader:
            raw_d = (row.get(sample_col) or "").strip()
            do = sp_to_donor.get(raw_d, raw_d)
            if do not in rna_donors:
                continue
            tmb_counts[do] = tmb_counts.get(do, 0) + 1
            ch = (row.get(chr_col) or "").replace("chr", "").strip()
            try:
                pos = int(float(row.get(pos_col, 0)))
                wgs_vars.append((do, ch, pos))
            except ValueError:
                continue

    matched_donors = sorted(tmb_counts.keys())
    print(f"[PCAWG-nonUS] Matched WGS + RNA-seq donors: {len(matched_donors)} ({len(wgs_vars):,} WGS SNVs)")

    pcawg_targets = [
        ("TERT", "ENSG00000164362", "5", 1295228, 30),
        ("PLEKHS1", "ENSG00000148735", "10", 115511590, 25),
        ("WDR74", "ENSG00000133316", "11", 62609080, 50),
        ("PAX5", "ENSG00000196092", "9", 37034000, 1500),
        ("BCL2", "ENSG00000171791", "18", 60985800, 1500),
        ("BCL6", "ENSG00000113916", "3", 187463000, 1500),
    ]
    ensg_lookup = {t[1]: t[0] for t in pcawg_targets}
    sym_lookup = {t[0]: t[0] for t in pcawg_targets}

    expr_map: dict[str, dict[str, float]] = {}
    with open_text_auto(expr_path) as f:
        cols = f.readline().rstrip("\n").split("\t")[1:]
        for line in f:
            idx = line.find("\t")
            if idx == -1:
                continue
            raw_g = line[:idx].split("|")[0].strip()
            base_g = raw_g.split(".")[0]
            gene = ensg_lookup.get(base_g) or sym_lookup.get(base_g)
            if not gene:
                continue
            vals = line.rstrip("\n").split("\t")[1:]
            dmap = {}
            for d, vs in zip(cols, vals):
                if d in tmb_counts:
                    try:
                        fv = float(vs)
                        if not math.isnan(fv):
                            dmap[d] = fv
                    except ValueError:
                        continue
            expr_map[gene] = dmap

    cnv_map: dict[str, dict[str, float]] = {t[0]: {} for t in pcawg_targets}
    if os.path.exists(cnv_path):
        with open_text_auto(cnv_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                samp = (row.get("sampleID") or "").strip()
                do = sp_to_donor.get(samp, samp if samp.startswith("DO") else "")
                if do not in tmb_counts:
                    continue
                ch = (row.get("chr") or "").replace("chr", "").strip()
                try:
                    st = int(float(row.get("start", 0)))
                    en = int(float(row.get("end", 0)))
                    cn = float(row.get("total_cn") or row.get("value") or 2.0)
                except ValueError:
                    continue
                for gene, _, t_ch, t_pos, _ in pcawg_targets:
                    if ch == t_ch and st <= t_pos <= en and do not in cnv_map[gene]:
                        cnv_map[gene][do] = cn

    pur_sorted = sorted(donor_purity[d] for d in matched_donors if d in donor_purity)
    med_pur = pur_sorted[len(pur_sorted) // 2] if pur_sorted else 0.70

    results = []
    for gene, _, chrom, center, tol in pcawg_targets:
        if gene not in expr_map:
            continue
        mut_set = {d for d, c, pos in wgs_vars if c == chrom and abs(pos - center) <= tol}
        mut_hists = {donor_hist.get(d, "") for d in mut_set if donor_hist.get(d, "")}
        X, y, wt_v, mut_v = [], [], [], []
        for d in matched_donors:
            if d not in expr_map[gene]:
                continue
            if mut_hists and donor_hist.get(d, "") not in mut_hists:
                continue
            yv = expr_map[gene][d]
            is_m = 1.0 if d in mut_set else 0.0
            cnv_v = cnv_map.get(gene, {}).get(d, 2.0)
            pur_v = donor_purity.get(d, med_pur)
            tmb_v = math.log10(tmb_counts.get(d, 0) + 1.0)
            X.append([1.0, is_m, cnv_v, pur_v, tmb_v])
            y.append(yv)
            (mut_v if is_m == 1.0 else wt_v).append(yv)

        if len(mut_v) < 3:
            continue

        beta, p_vals = fit_ols(X, y)
        n0, n1 = len(wt_v), len(mut_v)
        m0, m1 = sum(wt_v) / n0, sum(mut_v) / n1
        v0 = sum((v - m0) ** 2 for v in wt_v) / (n0 - 1)
        v1 = sum((v - m1) ** 2 for v in mut_v) / (n1 - 1)
        sp2 = ((n0 - 1) * v0 + (n1 - 1) * v1) / (n0 + n1 - 2)
        t_u = (m1 - m0) / math.sqrt(sp2 * (1.0 / n0 + 1.0 / n1))
        df_u = n0 + n1 - 2
        p_u = betainc(0.5 * df_u, 0.5, df_u / (df_u + t_u * t_u))

        results.append({
            "Dataset": f"PCAWG ({','.join(sorted(mut_hists)) or 'PanCan'})",
            "Gene": gene,
            "Locus_GRCh38": f"chr{chrom}:{center} (hg19)",
            "WT_Donors": n0,
            "Mut_Donors": n1,
            "Mean_WT_Expr": round(m0, 3),
            "Mean_Mut_Expr": round(m1, 3),
            "Unadj_Log2FC": round(m1 - m0, 3),
            "Unadj_P": p_u,
            "CovAdj_Beta_Mut": round(beta[1], 3),
            "CovAdj_P": p_vals[1],
        })
    return results


# ==============================================================================
# 4. Main Execution
# ==============================================================================

if __name__ == "__main__":
    all_rows = run_tcga_skcm_eqtl() + run_pcawg_nonus_eqtl()
    q_list = bh_fdr([r["CovAdj_P"] for r in all_rows])
    for r, q in zip(all_rows, q_list):
        r["Q_eQTL_BH"] = q

    out_path = os.path.join("results", "tables", "31_covariate_adjusted_eqtl_results.csv")
    if all_rows:
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"\nSaved genuine covariate-adjusted eQTL results -> {out_path}\n")
        for r in all_rows:
            print(
                f"{r['Dataset']:24s} | {r['Gene']:8s} | WT/MUT={r['WT_Donors']}/{r['Mut_Donors']} | "
                f"WT={r['Mean_WT_Expr']:.3f} vs MUT={r['Mean_Mut_Expr']:.3f} | "
                f"beta_mut={r['CovAdj_Beta_Mut']:+.3f} | p_unadj={r['Unadj_P']:.2e} | "
                f"p_cov={r['CovAdj_P']:.2e} | q_eQTL={r['Q_eQTL_BH']:.2e}"
            )
    else:
        print("No testable loci with >=3 mutated donors found.")