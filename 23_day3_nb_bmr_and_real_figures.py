import csv
import glob
import gzip
import math
import random
import statistics
from collections import Counter
from pathlib import Path


# ==============================================================================
# 1. Pure-Python Negative Binomial & Benjamini-Hochberg Helpers (Zero DLLs)
# ==============================================================================

LN10 = math.log(10.0)


def nb_log10_sf_exact(k: int, mu: float, alpha: float) -> float:
    """
    Exact upper-tail survival function P(X >= k) in -log10 space for
    Negative Binomial X ~ NB(mean=mu, dispersion=alpha), where Var(X) = mu + alpha * mu^2.
    """
    k = int(k)
    if k <= 0:
        return 0.0
    mu = max(float(mu), 1e-12)
    alpha = max(float(alpha), 1e-8)
    r = 1.0 / alpha
    p = r / (r + mu)
    q = mu / (r + mu)

    def logpmf(n: int) -> float:
        return (
            math.lgamma(n + r)
            - math.lgamma(r)
            - math.lgamma(n + 1)
            + r * math.log(p)
            + n * math.log(q)
        )

    if k > mu:
        s = 1.0
        t = 1.0
        for j in range(1, 250000):
            t *= (r + k + j - 1.0) / (k + j) * q
            s += t
            if t < 1e-16 * s:
                break
        logp = logpmf(k) + math.log(s)
    else:
        n = k - 1
        s = 1.0
        t = 1.0
        for j in range(1, n + 1):
            t *= (n - j + 1.0) / ((n - j + r) * q)
            s += t
            if t < 1e-16 * s:
                break
        cdf = min(math.exp(logpmf(n)) * s, 1.0)
        logp = math.log(max(1.0 - cdf, 1e-300))
    return float(min(max(-logp / LN10, 0.0), 300.0))


def poisson_log10_sf(k: int, mu: float) -> float:
    """Upper-tail Poisson P(X >= k) in -log10 space for direct baseline comparison."""
    k = int(k)
    if k <= 0:
        return 0.0
    mu = max(float(mu), 1e-12)
    if k > mu:
        log_pmf = -mu + k * math.log(mu) - math.lgamma(k + 1)
        s = 1.0
        t = 1.0
        for j in range(1, 200000):
            t *= mu / (k + j)
            s += t
            if t < 1e-16 * s:
                break
        logp = log_pmf + math.log(s)
    else:
        s = 1.0
        t = 1.0
        for j in range(1, k):
            t *= mu / j
            s += t
        cdf = min(math.exp(-mu) * s, 1.0)
        logp = math.log(max(1.0 - cdf, 1e-300))
    return float(min(max(-logp / LN10, 0.0), 300.0))


def bh_fdr_from_nl10(nl10_list: list[float], m_tests: int) -> list[float]:
    """Computes Benjamini-Hochberg q-values given -log10(p) values and total tests m_tests."""
    n = len(nl10_list)
    if n == 0:
        return []
    indexed = sorted(enumerate(nl10_list), key=lambda item: -item[1])
    q_out = [1.0] * n
    prev = 1.0
    for rank_idx in range(n - 1, -1, -1):
        orig_idx, nl = indexed[rank_idx]
        rank = rank_idx + 1
        p_val = 10.0 ** (-min(nl, 300.0))
        q = min(prev, p_val * float(m_tests) / float(rank))
        q_out[orig_idx] = min(max(q, 0.0), 1.0)
        prev = q
    return q_out


def open_text_auto(path: Path):
    with open(path, "rb") as raw:
        magic = raw.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


# ==============================================================================
# 2. Part A: Genome-Wide Zero-Included NB Dispersion Calibration (Melanoma WGS)
# ==============================================================================

def run_melanoma_nb_calibration(base: Path):
    snv_path = base / "results" / "pcawg_melanoma_wgs_snvs.csv.gz"
    clust_path = base / "results" / "tables" / "29_nb_clustered_noncoding_hotspots.csv"
    out_calib = base / "results" / "tables" / "32_melanoma_nb_bmr_calibration.csv"
    out_clusters = base / "results" / "tables" / "33_melanoma_recalibrated_nb_clusters.csv"

    print("=== PART A: Genome-Wide Zero-Included Negative Binomial BMR Calibration ===")
    valid_chroms = {str(i) for i in range(1, 23)} | {"X", "Y"}
    window_size = 50

    donors: set[str] = set()
    chrom_max: dict[str, int] = {}
    bin_counts: Counter = Counter()
    n_snvs = 0

    with gzip.open(snv_path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ch = (row.get("Chromosome") or "").replace("chr", "").strip()
            if ch not in valid_chroms:
                continue
            pos_raw = (row.get("Start_Position") or "").strip()
            if not pos_raw or pos_raw.lower() == "nan":
                continue
            try:
                pos = int(float(pos_raw))
            except ValueError:
                continue
            pat = (row.get("patient") or "").strip()
            if pat:
                donors.add(pat)
            n_snvs += 1
            if pos > chrom_max.get(ch, 0):
                chrom_max[ch] = pos
            bin_counts[(ch, pos // window_size)] += 1

    n_donors = len(donors)
    total_windows = sum(math.ceil(mx / window_size) for mx in chrom_max.values())
    n_nonempty = len(bin_counts)
    n_zero = total_windows - n_nonempty

    mean_all = float(n_snvs) / float(total_windows)
    sum_sq_all = sum(float(c) * float(c) for c in bin_counts.values())
    var_all = (sum_sq_all - total_windows * (mean_all ** 2)) / float(total_windows - 1)
    alpha_zero_inc = max((var_all - mean_all) / (mean_all ** 2), 1e-6)
    shape_zero_inc = 1.0 / alpha_zero_inc

    mean_trunc = float(n_snvs) / float(n_nonempty)
    var_trunc = (sum_sq_all - n_nonempty * (mean_trunc ** 2)) / float(n_nonempty - 1)
    alpha_trunc = max((var_trunc - mean_trunc) / (mean_trunc ** 2), 1e-6)

    print(f"Cohort: PCAWG Melanoma WGS (N = {n_donors} donors, {n_snvs:,} SNVs)")
    print(f"Total 50-bp genome-wide windows (M_total): {total_windows:,} (Zero-count: {n_zero:,}, Non-zero: {n_nonempty:,})")
    print(f"Zero-Included Mean (mu_0):   {mean_all:.5f} SNVs / 50 bp")
    print(f"Zero-Included Variance:      {var_all:.5f} (Fano factor = {var_all / mean_all:.2f})")
    print(f"Zero-Included NB alpha_hat:  {alpha_zero_inc:.4f} (shape r = {shape_zero_inc:.4f})")
    print(f"Naive Zero-Truncated alpha:  {alpha_trunc:.4f} (biased if Obs=0 windows omitted)")

    calib_rows = [
        {
            "Metric": "Zero_Included_GenomeWide_NB",
            "Donors": n_donors,
            "Total_SNVs": n_snvs,
            "Window_BP": window_size,
            "Total_Windows_M": total_windows,
            "Zero_Windows": n_zero,
            "NonZero_Windows": n_nonempty,
            "Mean_SNVs_Per_Window": round(mean_all, 6),
            "Variance_Per_Window": round(var_all, 6),
            "Fano_Factor_Var_Over_Mean": round(var_all / mean_all, 4),
            "NB_Dispersion_Alpha": round(alpha_zero_inc, 4),
            "NB_Shape_r": round(shape_zero_inc, 4),
        },
        {
            "Metric": "Naive_ZeroTruncated_ObsGe1",
            "Donors": n_donors,
            "Total_SNVs": n_snvs,
            "Window_BP": window_size,
            "Total_Windows_M": n_nonempty,
            "Zero_Windows": 0,
            "NonZero_Windows": n_nonempty,
            "Mean_SNVs_Per_Window": round(mean_trunc, 6),
            "Variance_Per_Window": round(var_trunc, 6),
            "Fano_Factor_Var_Over_Mean": round(var_trunc / mean_trunc, 4),
            "NB_Dispersion_Alpha": round(alpha_trunc, 4),
            "NB_Shape_r": round(1.0 / alpha_trunc, 4),
        },
    ]
    with open(out_calib, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(calib_rows[0].keys()))
        writer.writeheader()
        writer.writerows(calib_rows)
    print(f"Saved calibration summary -> {out_calib}")

    if clust_path.exists():
        with open(clust_path, "r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        nl10_nb = []
        nl10_pois = []
        for r in rows:
            obs = int(float(r["Observed_Mutations"]))
            exp_loc = float(r["Local_Expected"])
            nl10_nb.append(nb_log10_sf_exact(obs, exp_loc, alpha_zero_inc))
            nl10_pois.append(poisson_log10_sf(obs, exp_loc))

        q_nb = bh_fdr_from_nl10(nl10_nb, total_windows)
        q_pois = bh_fdr_from_nl10(nl10_pois, total_windows)

        for r, pnb, qn, qp in zip(rows, nl10_nb, q_nb, q_pois):
            r["NB_Alpha_ZeroInc"] = round(alpha_zero_inc, 4)
            r["NegLog10_P_NB_ZeroInc"] = round(pnb, 3)
            r["Q_NB_ZeroInc"] = f"{qn:.3e}"
            r["NegLog10_Q_NB_ZeroInc"] = round(-math.log10(max(qn, 1e-300)), 3)
            r["Q_Poisson_Uncalibrated"] = f"{qp:.3e}"
            r["_q_sort"] = qn

        rows.sort(key=lambda x: x["_q_sort"])
        for r in rows:
            del r["_q_sort"]

        with open(out_clusters, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved recalibrated cluster table -> {out_clusters}\n")
        print("Top 10 Melanoma Noncoding Clusters (Zero-Included NB BMR):")
        print(f"{'Chr':4s} | {'Peak_Center':11s} | {'Target_Gene':14s} | {'Obs':4s} | {'Pat':3s} | {'LocExp':7s} | {'Q_NB_ZeroInc':12s} | {'Q_Poisson':12s}")
        print("-" * 82)
        for r in rows[:10]:
            print(
                f"{r['Chromosome']:4s} | {r['Peak_Center']:11s} | {r['Target_Gene'][:14]:14s} | "
                f"{int(float(r['Observed_Mutations'])):4d} | {int(float(r['Mutated_Patients'])):3d} | "
                f"{float(r['Local_Expected']):7.3f} | {r['Q_NB_ZeroInc']:12s} | {r['Q_Poisson_Uncalibrated']:12s}"
            )


# ==============================================================================
# 3. Part B: Regenerate Publication eQTL Figure from Genuine PCAWG Data
# ==============================================================================

def _quantile(sorted_vals: list[float], q: float) -> float:
    if not sorted_vals:
        return 0.0
    pos = q * (len(sorted_vals) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_vals[lo]
    return sorted_vals[lo] + (pos - lo) * (sorted_vals[hi] - sorted_vals[lo])


def _write_vector_svg_eqtl(
    svg_path: Path,
    panels: list[tuple[str, str, str, int, list[float], list[float], float, float, float]],
):
    W, H = 1100, 860
    rng = random.Random(42)
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">',
        f'<text x="{W/2}" y="34" text-anchor="middle" font-size="18" font-weight="bold" fill="#111111">Covariate-Adjusted Cis-Regulatory eQTL Validation in PCAWG WGS + RNA-seq (N = 433 Donors)</text>',
        f'<text x="{W/2}" y="56" text-anchor="middle" font-size="13" fill="#444444">OLS Model: log2(FPKM-UQ + 1) ~ Mut_Status + Total_CN + Tumor_Purity + log10(TMB + 1)</text>',
    ]

    coords = [(75, 95), (585, 95), (75, 480), (585, 480)]
    pw, ph = 445, 315

    for idx, (gene, region_lbl, chrom, pos, wt_v, mut_v, beta_m, p_cov, q_val) in enumerate(panels):
        ox, oy = coords[idx]
        is_sig = gene != "WDR74"
        all_v = wt_v + mut_v
        ymin = 0.0
        ymax = max(max(all_v) * 1.18, 1.0)

        def y_px(val: float) -> float:
            return oy + ph - (val - ymin) / (ymax - ymin) * ph

        lines.append(f'<text x="{ox + pw/2}" y="{oy - 18}" text-anchor="middle" font-size="14" font-weight="bold" fill="#111111">{gene} ({region_lbl})</text>')
        lines.append(f'<text x="{ox + pw/2}" y="{oy - 4}" text-anchor="middle" font-size="11.5" fill="#555555">chr{chrom}:{pos:,} (hg19)</text>')
        lines.append(f'<line x1="{ox}" y1="{oy}" x2="{ox}" y2="{oy+ph}" stroke="#222222" stroke-width="1.3"/>')
        lines.append(f'<line x1="{ox}" y1="{oy+ph}" x2="{ox+pw}" y2="{oy+ph}" stroke="#222222" stroke-width="1.3"/>')
        lines.append(f'<text x="{ox - 42}" y="{oy + ph/2}" text-anchor="middle" font-size="12" fill="#222222" transform="rotate(-90 {ox - 42},{oy + ph/2})">log2(FPKM-UQ + 1)</text>')

        n_ticks = 5
        for t in range(n_ticks + 1):
            tv = ymin + t * (ymax - ymin) / n_ticks
            yp = y_px(tv)
            lines.append(f'<line x1="{ox-5}" y1="{yp:.1f}" x2="{ox}" y2="{yp:.1f}" stroke="#222222" stroke-width="1.1"/>')
            lines.append(f'<text x="{ox-9}" y="{yp+4:.1f}" text-anchor="end" font-size="10.5" fill="#333333">{tv:.1f}</text>')

        groups = [
            (ox + pw * 0.30, sorted(wt_v), f"Wild-Type (n={len(wt_v)})", "#4C72B0"),
            (ox + pw * 0.72, sorted(mut_v), f"Mutated (n={len(mut_v)})", "#DD5143" if is_sig else "#7F8C8D"),
        ]
        for gx, svals, glbl, col in groups:
            q25, q50, q75 = _quantile(svals, 0.25), _quantile(svals, 0.50), _quantile(svals, 0.75)
            iqr = q75 - q25
            w_lo = min([v for v in svals if v >= q25 - 1.5 * iqr] or [svals[0]])
            w_hi = max([v for v in svals if v <= q75 + 1.5 * iqr] or [svals[-1]])
            bw = 56
            lines.append(f'<line x1="{gx}" y1="{y_px(w_lo):.1f}" x2="{gx}" y2="{y_px(w_hi):.1f}" stroke="#222222" stroke-width="1.2"/>')
            lines.append(f'<line x1="{gx-14}" y1="{y_px(w_lo):.1f}" x2="{gx+14}" y2="{y_px(w_lo):.1f}" stroke="#222222" stroke-width="1.2"/>')
            lines.append(f'<line x1="{gx-14}" y1="{y_px(w_hi):.1f}" x2="{gx+14}" y2="{y_px(w_hi):.1f}" stroke="#222222" stroke-width="1.2"/>')
            y_top, y_bot = y_px(q75), y_px(q25)
            lines.append(f'<rect x="{gx - bw/2}" y="{y_top:.1f}" width="{bw}" height="{max(y_bot - y_top, 1.0):.1f}" fill="{col}" fill-opacity="0.30" stroke="#222222" stroke-width="1.2"/>')
            lines.append(f'<line x1="{gx - bw/2}" y1="{y_px(q50):.1f}" x2="{gx + bw/2}" y2="{y_px(q50):.1f}" stroke="#111111" stroke-width="2.2"/>')
            for v in svals:
                jx = gx + rng.uniform(-14.0, 14.0)
                lines.append(f'<circle cx="{jx:.1f}" cy="{y_px(v):.1f}" r="2.7" fill="{col}" fill-opacity="0.65"/>')
            lines.append(f'<text x="{gx}" y="{oy + ph + 20}" text-anchor="middle" font-size="11.5" font-weight="bold" fill="#222222">{glbl}</text>')

        bg_c = "#FFF7F5" if is_sig else "#F4F6F8"
        bd_c = "#DD5143" if is_sig else "#AAB2BD"
        lines.append(f'<rect x="{ox + 12}" y="{oy + 10}" width="248" height="44" rx="5" fill="{bg_c}" stroke="{bd_c}" stroke-width="1.1"/>')
        lines.append(f'<text x="{ox + 20}" y="{oy + 27}" font-size="11" font-weight="bold" fill="#111111">Cov-Adj beta_mut = {beta_m:+.3f}</text>')
        lines.append(f'<text x="{ox + 20}" y="{oy + 44}" font-size="10.5" fill="#333333">p_cov = {p_cov:.2e} | q_eQTL = {q_val:.2e}</text>')

    lines.append("</svg>")
    svg_path.write_text("\n".join(lines), encoding="utf-8")


def regenerate_real_eqtl_figure(base: Path):
    print("\n=== PART B: Regenerating Real PCAWG Covariate-Adjusted eQTL Figure ===")
    eqtl_csv = base / "results" / "tables" / "31_covariate_adjusted_eqtl_results.csv"
    expr_path = base / "results" / "external_pcawg" / "pcawg_rnaseq_fpkm_uq.tsv.gz"
    sp_hist = base / "pcawg_specimen_histology_August2016_v9"
    do_hist = base / "pcawg_specimen_histology_August2016_v9_donor"
    maf_hits = [
        Path(p) for p in glob.glob(str(base / "**" / "*whitelist*.nonUS*"), recursive=True)
        if "\\venv\\" not in p and "/venv/" not in p
    ]

    if not (eqtl_csv.exists() and expr_path.exists() and maf_hits):
        print("[Skip Part B] Required cached PCAWG files not found.")
        return

    with open(eqtl_csv, "r", encoding="utf-8") as f:
        stats_by_gene = {r["Gene"]: r for r in csv.DictReader(f)}

    sp_to_donor = {}
    donor_hist = {}
    with open_text_auto(sp_hist) as f_sp, open_text_auto(do_hist) as f_do:
        for rsp, rdo in zip(csv.DictReader(f_sp, delimiter="\t"), csv.DictReader(f_do, delimiter="\t")):
            sp = (rsp.get("icgc_specimen_id") or "").strip()
            do = (rdo.get("icgc_specimen_id") or "").strip()
            hi = (rdo.get("histology_abbreviation") or "").strip()
            if sp and do:
                sp_to_donor[sp] = do
                sp_to_donor[do] = do
            if do and hi:
                donor_hist[do] = hi

    with open_text_auto(expr_path) as f:
        rna_donors = set(f.readline().rstrip("\n").split("\t")[1:])

    targets = [
        ("BCL6", "ENSG00000113916", "3", 187463000, 1500, "5'-UTR / Intron 1 SE"),
        ("BCL2", "ENSG00000171791", "18", 60985800, 1500, "Promoter / 5'-UTR"),
        ("PAX5", "ENSG00000196092", "9", 37034000, 1500, "Enhancer / TSS-proximal"),
        ("WDR74", "ENSG00000133316", "11", 62609080, 50, "Promoter (Negative Control)"),
    ]
    ensg_lookup = {t[1]: t[0] for t in targets}
    sym_lookup = {t[0]: t[0] for t in targets}

    mut_donors = {t[0]: set() for t in targets}
    with open_text_auto(maf_hits[0]) as f:
        reader = csv.DictReader(f, delimiter="\t")
        scol = next(c for c in reader.fieldnames if any(k in c.lower() for k in ("sample", "donor", "patient")))
        ccol = next(c for c in reader.fieldnames if c.lower() in ("chr", "chrom", "chromosome", "#chrom"))
        pcol = next(c for c in reader.fieldnames if "start" in c.lower())
        for row in reader:
            do = sp_to_donor.get((row.get(scol) or "").strip(), (row.get(scol) or "").strip())
            if do not in rna_donors:
                continue
            ch = (row.get(ccol) or "").replace("chr", "").strip()
            try:
                pos = int(float(row.get(pcol, 0)))
            except ValueError:
                continue
            for gene, _, t_ch, t_pos, tol, _ in targets:
                if ch == t_ch and abs(pos - t_pos) <= tol:
                    mut_donors[gene].add(do)

    expr_map = {}
    with open_text_auto(expr_path) as f:
        cols = f.readline().rstrip("\n").split("\t")[1:]
        for line in f:
            idx = line.find("\t")
            if idx == -1:
                continue
            raw_g = line[:idx].split("|")[0].strip().split(".")[0]
            gene = ensg_lookup.get(raw_g) or sym_lookup.get(raw_g)
            if not gene:
                continue
            vals = line.rstrip("\n").split("\t")[1:]
            dmap = {}
            for d, vs in zip(cols, vals):
                if d in rna_donors:
                    try:
                        fv = float(vs)
                        if not math.isnan(fv):
                            dmap[d] = fv
                    except ValueError:
                        continue
            expr_map[gene] = dmap

    panels = []
    for gene, _, chrom, pos, _, region_lbl in targets:
        mset = mut_donors[gene]
        m_hists = {donor_hist.get(d, "") for d in mset if donor_hist.get(d, "")}
        wt_vals, mut_vals = [], []
        for d, val in expr_map.get(gene, {}).items():
            if m_hists and donor_hist.get(d, "") not in m_hists:
                continue
            if d in mset:
                mut_vals.append(val)
            else:
                wt_vals.append(val)
        st = stats_by_gene[gene]
        panels.append((
            gene, region_lbl, chrom, pos, wt_vals, mut_vals,
            float(st["CovAdj_Beta_Mut"]), float(st["CovAdj_P"]), float(st["Q_eQTL_BH"])
        ))

    fig_dir = base / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    svg_out = fig_dir / "fig17_eqtl_expression_correlation.svg"
    _write_vector_svg_eqtl(svg_out, panels)
    print(f"Saved publication vector SVG eQTL figure -> {svg_out}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plt.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        })
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0), dpi=300)
        axes = axes.flatten()
        rng = random.Random(42)

        for ax, (gene, region_lbl, chrom, pos, wt_vals, mut_vals, beta_m, p_cov, q_val) in zip(axes, panels):
            bp = ax.boxplot(
                [wt_vals, mut_vals],
                positions=[1, 2],
                widths=0.45,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "#111111", "linewidth": 2.0},
            )
            is_sig = gene != "WDR74"
            colors = ["#4C72B0", "#DD5143" if is_sig else "#8C969E"]
            for patch, col in zip(bp["boxes"], colors):
                patch.set_facecolor(col)
                patch.set_alpha(0.35)
            ax.scatter([1.0 + rng.uniform(-0.08, 0.08) for _ in wt_vals], wt_vals, s=14, color="#4C72B0", alpha=0.55, zorder=3)
            ax.scatter([2.0 + rng.uniform(-0.08, 0.08) for _ in mut_vals], mut_vals, s=22, color="#DD5143" if is_sig else "#555555", alpha=0.80, zorder=4)
            ax.set_xticks([1, 2])
            ax.set_xticklabels([f"Wild-Type\n(n = {len(wt_vals)})", f"Mutated\n(n = {len(mut_vals)})"])
            ax.set_ylabel("Expression [log2(FPKM-UQ + 1)]")
            ax.set_title(f"{gene} ({region_lbl})\nchr{chrom}:{pos:,} (hg19)", fontweight="bold")
            ax.text(
                0.04, 0.95,
                f"Cov-Adj beta_mut = {beta_m:+.3f}\np_cov = {p_cov:.2e} | q_eQTL = {q_val:.2e}",
                transform=ax.transAxes,
                fontsize=8.8,
                verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.35", facecolor="#FFF8F6" if is_sig else "#F4F6F8", edgecolor="#DD5143" if is_sig else "#B0B8C0"),
            )
        fig.suptitle(
            "Covariate-Adjusted Cis-Regulatory eQTL Validation in PCAWG WGS + RNA-seq (N = 433 Donors)\n"
            "OLS Model: log2(FPKM-UQ + 1) ~ Mut_Status + Total_CN + Tumor_Purity + log10(TMB + 1)",
            fontsize=11.5, fontweight="bold", y=0.98
        )
        plt.tight_layout(rect=[0, 0, 1, 0.94])
        png_out = fig_dir / "fig17_eqtl_expression_correlation.png"
        pdf_out = fig_dir / "fig17_eqtl_expression_correlation.pdf"
        fig.savefig(png_out, dpi=300)
        fig.savefig(pdf_out)
        plt.close(fig)
        print(f"Saved high-res PNG & PDF eQTL figure    -> {png_out}")
    except Exception as exc:
        print(f"[Note] Matplotlib raster export skipped ({exc}); vector SVG saved cleanly.")


# ==============================================================================
# 4. Part C: Pure-Python Audit of ML Benchmark Matrix
# ==============================================================================

def audit_ml_benchmark_matrix(base: Path):
    print("\n=== PART C: Auditing ML Hotspot Prioritization Matrix ===")
    ml_path = base / "results" / "tables" / "ml_hotspot_training_matrix.csv"
    if not ml_path.exists():
        print("[Skip Part C] ml_hotspot_training_matrix.csv not found.")
        return
    with open(ml_path, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    print(f"Total hotspot loci in matrix: {len(rows)}")
    by_cancer = Counter(r["Cancer_Type"] for r in rows)
    for k, v in by_cancer.items():
        print(f"  - {k}: {v} loci")

    num_cols = ["Mutated_Donors", "Fold_Enrichment", "Absolute_Distance_to_TSS", "Is_Promoter", "Is_Enhancer", "Has_TF_Impact"]
    print("\nFeature Summary (Mean +/- SD | Min - Max):")
    for col in num_cols:
        vals = [float(r[col]) for r in rows if r.get(col) not in ("", None)]
        if vals:
            mean_v = statistics.mean(vals)
            sd_v = statistics.stdev(vals) if len(vals) > 1 else 0.0
            print(f"  {col:26s}: {mean_v:10.3f} +/- {sd_v:9.3f}  |  [{min(vals):.3f} .. {max(vals):.3f}]")


if __name__ == "__main__":
    BASE_DIR = Path(r"D:\noncoding mutation finding surver")
    run_melanoma_nb_calibration(BASE_DIR)
    regenerate_real_eqtl_figure(BASE_DIR)
    audit_ml_benchmark_matrix(BASE_DIR)