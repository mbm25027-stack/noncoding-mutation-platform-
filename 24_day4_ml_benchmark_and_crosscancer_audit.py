import csv
import gzip
import math
from collections import Counter
from pathlib import Path


# ==============================================================================
# 1. Pure-Python NB, FDR, ROC-AUC, PR-AUC & LOOCV Classifiers (Zero DLLs)
# ==============================================================================

LN10 = math.log(10.0)


def nb_log10_sf_exact(k: int, mu: float, alpha: float) -> float:
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


def compute_roc_pr_curves(y_true: list[int], y_score: list[float]):
    """Exact ROC-AUC and Average Precision (PR-AUC) without external libraries."""
    n_pos = sum(1 for y in y_true if y == 1)
    n_neg = len(y_true) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5, 0.0, [(0.0, 0.0), (1.0, 1.0)], [(0.0, 1.0), (1.0, 0.0)]

    paired = sorted(zip(y_score, y_true), key=lambda x: -x[0])
    tp = 0
    fp = 0
    roc_pts = [(0.0, 0.0)]
    pr_pts = [(0.0, 1.0)]
    ap = 0.0

    i = 0
    n = len(paired)
    while i < n:
        j = i
        while j < n and paired[j][0] == paired[i][0]:
            j += 1
        batch_pos = sum(paired[k][1] for k in range(i, j))
        batch_neg = (j - i) - batch_pos
        prev_rec = tp / n_pos
        tp += batch_pos
        fp += batch_neg
        rec = tp / n_pos
        fpr = fp / n_neg
        prec = tp / (tp + fp)
        roc_pts.append((fpr, rec))
        pr_pts.append((rec, prec))
        ap += (rec - prev_rec) * prec
        i = j

    roc_auc = 0.0
    for idx in range(1, len(roc_pts)):
        x0, y0 = roc_pts[idx - 1]
        x1, y1 = roc_pts[idx]
        roc_auc += (x1 - x0) * 0.5 * (y0 + y1)

    return roc_auc, ap, roc_pts, pr_pts


def fit_predict_logistic_l2(
    X_train: list[list[float]],
    y_train: list[int],
    x_test: list[float],
    l2: float = 0.5,
    lr: float = 0.15,
    epochs: int = 400,
) -> float:
    """In-fold standardized L2 Logistic Regression (zero test-set leakage)."""
    n = len(X_train)
    d = len(x_test)
    means = [sum(X_train[i][j] for i in range(n)) / n for j in range(d)]
    stds = []
    for j in range(d):
        var = sum((X_train[i][j] - means[j]) ** 2 for i in range(n)) / max(n - 1, 1)
        stds.append(math.sqrt(var) if var > 1e-12 else 1.0)

    Xs = [[(X_train[i][j] - means[j]) / stds[j] for j in range(d)] for i in range(n)]
    xt = [(x_test[j] - means[j]) / stds[j] for j in range(d)]

    w = [0.0] * d
    b = 0.0
    for _ in range(epochs):
        gw = [0.0] * d
        gb = 0.0
        for i in range(n):
            z = b + sum(w[j] * Xs[i][j] for j in range(d))
            z = max(min(z, 30.0), -30.0)
            p = 1.0 / (1.0 + math.exp(-z))
            err = p - y_train[i]
            gb += err
            for j in range(d):
                gw[j] += err * Xs[i][j]
        b -= lr * (gb / n)
        for j in range(d):
            w[j] -= lr * ((gw[j] / n) + l2 * w[j])

    zt = max(min(b + sum(w[j] * xt[j] for j in range(d)), 30.0), -30.0)
    return 1.0 / (1.0 + math.exp(-zt))


def fit_predict_gbdt_stumps(
    X_train: list[list[float]],
    y_train: list[int],
    x_test: list[float],
    n_estimators: int = 18,
    lr: float = 0.18,
) -> float:
    """Pure-Python Gradient Boosted Decision Stumps (log-loss)."""
    n = len(X_train)
    d = len(x_test)
    pos_rate = min(max(sum(y_train) / n, 1e-4), 1.0 - 1e-4)
    f0 = math.log(pos_rate / (1.0 - pos_rate))
    F = [f0] * n
    pred_test = f0

    for _ in range(n_estimators):
        res = []
        hess = []
        for i in range(n):
            z = max(min(F[i], 30.0), -30.0)
            p = 1.0 / (1.0 + math.exp(-z))
            res.append(y_train[i] - p)
            hess.append(max(p * (1.0 - p), 1e-6))

        best_gain = -1.0
        best_stump = None
        for j in range(d):
            vals = sorted(set(X_train[i][j] for i in range(n)))
            if len(vals) <= 1:
                continue
            thresholds = [0.5 * (vals[k] + vals[k + 1]) for k in range(len(vals) - 1)]
            for thr in thresholds:
                g_l = h_l = g_r = h_r = 0.0
                for i in range(n):
                    if X_train[i][j] <= thr:
                        g_l += res[i]
                        h_l += hess[i]
                    else:
                        g_r += res[i]
                        h_r += hess[i]
                lam = 1.5
                gain = 0.5 * ((g_l * g_l) / (h_l + lam) + (g_r * g_r) / (h_r + lam))
                if gain > best_gain:
                    best_gain = gain
                    best_stump = (j, thr, g_l / (h_l + lam), g_r / (h_r + lam))

        if best_stump is None:
            break
        j_b, thr_b, val_l, val_r = best_stump
        for i in range(n):
            F[i] += lr * (val_l if X_train[i][j_b] <= thr_b else val_r)
        pred_test += lr * (val_l if x_test[j_b] <= thr_b else val_r)

    zt = max(min(pred_test, 30.0), -30.0)
    return 1.0 / (1.0 + math.exp(-zt))


# ==============================================================================
# 2. Part A: Cross-Cancer PCAWG Breast WGS (N = 120) NB BMR Calibration
# ==============================================================================

def run_breast_nb_calibration(base: Path):
    snv_path = base / "results" / "pcawg_breast_wgs_snvs.csv.gz"
    cand_path = base / "results" / "tables" / "breast_noncoding_candidates_annotated.csv"
    out_calib = base / "results" / "tables" / "34_breast_nb_bmr_calibration.csv"
    out_cands = base / "results" / "tables" / "35_breast_recalibrated_candidates.csv"

    print("=== PART A: Cross-Cancer PCAWG Breast WGS (N = 120) NB BMR Calibration ===")
    if not snv_path.exists():
        print("[Skip Part A] pcawg_breast_wgs_snvs.csv.gz not found.")
        return

    valid_chroms = {str(i) for i in range(1, 23)} | {"X", "Y"}
    window_size = 50
    donors: set[str] = set()
    chrom_max: dict[str, int] = {}
    bin_counts: Counter = Counter()
    n_snvs = 0

    with gzip.open(snv_path, "rt", encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f):
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

    print(f"Cohort: PCAWG Breast WGS (N = {n_donors} donors, {n_snvs:,} SNVs)")
    print(f"Total 50-bp windows (M_total): {total_windows:,} (Zero: {n_zero:,}, Non-zero: {n_nonempty:,})")
    print(f"Zero-Included Mean (mu_0):     {mean_all:.5f} SNVs / 50 bp")
    print(f"Zero-Included Variance:        {var_all:.5f} (Fano factor = {var_all / mean_all:.2f})")
    print(f"Zero-Included NB alpha_hat:    {alpha_zero_inc:.4f} (shape r = {shape_zero_inc:.4f})")

    with open(out_calib, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Cohort", "Donors", "Total_SNVs", "Window_BP", "Total_Windows_M",
                "Zero_Windows", "NonZero_Windows", "Mean_SNVs_Per_Window",
                "Variance_Per_Window", "Fano_Factor", "NB_Dispersion_Alpha", "NB_Shape_r"
            ],
        )
        writer.writeheader()
        writer.writerow({
            "Cohort": "PCAWG_Breast_WGS",
            "Donors": n_donors,
            "Total_SNVs": n_snvs,
            "Window_BP": window_size,
            "Total_Windows_M": total_windows,
            "Zero_Windows": n_zero,
            "NonZero_Windows": n_nonempty,
            "Mean_SNVs_Per_Window": round(mean_all, 6),
            "Variance_Per_Window": round(var_all, 6),
            "Fano_Factor": round(var_all / mean_all, 4),
            "NB_Dispersion_Alpha": round(alpha_zero_inc, 4),
            "NB_Shape_r": round(shape_zero_inc, 4),
        })
    print(f"Saved Breast NB calibration -> {out_calib}")

    if cand_path.exists():
        with open(cand_path, "r", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        nl10_nb = []
        nl10_pois = []
        for r in rows:
            obs = int(float(r["Mutated_Donors"]))
            fe = max(float(r.get("Fold_Enrichment", 1.0)), 1e-3)
            exp_loc = max(obs / fe, mean_all)
            r["_exp_loc"] = round(exp_loc, 4)
            nl10_nb.append(nb_log10_sf_exact(obs, exp_loc, alpha_zero_inc))
            nl10_pois.append(poisson_log10_sf(obs, exp_loc))

        q_nb = bh_fdr_from_nl10(nl10_nb, total_windows)
        q_pois = bh_fdr_from_nl10(nl10_pois, total_windows)
        for r, qn, qp in zip(rows, q_nb, q_pois):
            r["Est_Local_Expected"] = r.pop("_exp_loc")
            r["Q_NB_ZeroInc"] = f"{qn:.3e}"
            r["Q_Poisson_Uncalibrated"] = f"{qp:.3e}"
            r["_q"] = qn
        rows.sort(key=lambda x: x["_q"])
        for r in rows:
            del r["_q"]

        with open(out_cands, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved recalibrated Breast candidates -> {out_cands}\n")
        print("Top PCAWG Breast Noncoding Candidates (Zero-Included NB BMR):")
        print(f"{'Chr':4s} | {'Position':11s} | {'Target_Gene':14s} | {'Element':16s} | {'Donors':6s} | {'FE':7s} | {'Q_NB_ZeroInc':12s}")
        print("-" * 82)
        for r in rows[:10]:
            print(
                f"{str(r['Chromosome']).replace('chr',''):4s} | {r['Start_Position']:11s} | "
                f"{r['Target_Gene'][:14]:14s} | {r['Element_Class'][:16]:16s} | "
                f"{int(float(r['Mutated_Donors'])):6d} | {float(r['Fold_Enrichment']):7.1f} | {r['Q_NB_ZeroInc']:12s}"
            )


# ==============================================================================
# 3. Part B: Leakage-Free LOOCV ML Benchmark (PR-AUC & ROC-AUC)
# ==============================================================================

def _write_svg_pr_roc(
    svg_path: Path,
    curves: dict[str, tuple[float, float, list[tuple[float, float]], list[tuple[float, float]], str]],
):
    """Generates publication-grade 2-panel (ROC + PR) vector SVG figure."""
    W, H = 1080, 480
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">',
        f'<text x="{W/2}" y="30" text-anchor="middle" font-size="16.5" font-weight="bold" fill="#111111">Leakage-Free Leave-One-Out Cross-Validation (LOOCV) Prioritization Benchmark (N = 25 Loci)</text>',
    ]

    panels = [
        ("Receiver Operating Characteristic (ROC)", "False Positive Rate (1 - Specificity)", "True Positive Rate (Sensitivity)", 75, 65, True),
        ("Precision-Recall Curve (PR)", "Recall (Sensitivity)", "Precision (Positive Predictive Value)", 585, 65, False),
    ]
    pw, ph = 425, 345

    for title, xlab, ylab, ox, oy, is_roc in panels:
        lines.append(f'<text x="{ox + pw/2}" y="{oy - 12}" text-anchor="middle" font-size="13.5" font-weight="bold" fill="#111111">{title}</text>')
        lines.append(f'<line x1="{ox}" y1="{oy}" x2="{ox}" y2="{oy+ph}" stroke="#222222" stroke-width="1.3"/>')
        lines.append(f'<line x1="{ox}" y1="{oy+ph}" x2="{ox+pw}" y2="{oy+ph}" stroke="#222222" stroke-width="1.3"/>')
        lines.append(f'<text x="{ox + pw/2}" y="{oy + ph + 38}" text-anchor="middle" font-size="11.5" fill="#222222">{xlab}</text>')
        lines.append(f'<text x="{ox - 42}" y="{oy + ph/2}" text-anchor="middle" font-size="11.5" fill="#222222" transform="rotate(-90 {ox - 42},{oy + ph/2})">{ylab}</text>')

        for t in range(6):
            v = t * 0.2
            xp = ox + v * pw
            yp = oy + ph - v * ph
            lines.append(f'<line x1="{xp:.1f}" y1="{oy+ph}" x2="{xp:.1f}" y2="{oy+ph+4}" stroke="#222222" stroke-width="1"/>')
            lines.append(f'<text x="{xp:.1f}" y="{oy+ph+17}" text-anchor="middle" font-size="10" fill="#333333">{v:.1f}</text>')
            lines.append(f'<line x1="{ox-4}" y1="{yp:.1f}" x2="{ox}" y2="{yp:.1f}" stroke="#222222" stroke-width="1"/>')
            lines.append(f'<text x="{ox-8}" y="{yp+3.5:.1f}" text-anchor="end" font-size="10" fill="#333333">{v:.1f}</text>')

        if is_roc:
            lines.append(f'<line x1="{ox}" y1="{oy+ph}" x2="{ox+pw}" y2="{oy}" stroke="#999999" stroke-width="1.1" stroke-dasharray="5,4"/>')

        for idx_m, (m_name, (roc_auc, pr_auc, roc_pts, pr_pts, col)) in enumerate(curves.items()):
            pts = roc_pts if is_roc else pr_pts
            poly = " ".join(f"{ox + x * pw:.1f},{oy + ph - y * ph:.1f}" for x, y in pts)
            dash = ' stroke-dasharray="6,3"' if "Baseline" in m_name else ""
            lines.append(f'<polyline fill="none" stroke="{col}" stroke-width="2.4"{dash} points="{poly}"/>')

            ly = oy + ph - 78 + idx_m * 20
            lx = ox + 115 if is_roc else ox + 15
            metric_txt = f"ROC-AUC = {roc_auc:.3f}" if is_roc else f"PR-AUC = {pr_auc:.3f}"
            lines.append(f'<line x1="{lx}" y1="{ly-4}" x2="{lx+22}" y2="{ly-4}" stroke="{col}" stroke-width="2.4"{dash}/>')
            lines.append(f'<text x="{lx+28}" y="{ly}" font-size="10.5" fill="#111111">{m_name} ({metric_txt})</text>')

    lines.append("</svg>")
    svg_path.write_text("\n".join(lines), encoding="utf-8")


def run_ml_loocv_benchmark(base: Path):
    print("\n=== PART B: Leakage-Free LOOCV ML Benchmark (PR-AUC & ROC-AUC) ===")
    ml_path = base / "results" / "tables" / "ml_hotspot_training_matrix.csv"
    mel_path = base / "results" / "tables" / "melanoma_noncoding_candidates_annotated.csv"
    brc_path = base / "results" / "tables" / "breast_noncoding_candidates_annotated.csv"

    if not ml_path.exists():
        print("[Skip Part B] ml_hotspot_training_matrix.csv not found.")
        return

    locus_meta: dict[str, tuple[str, str]] = {}
    for p in (mel_path, brc_path):
        if not p.exists():
            continue
        with open(p, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                ch = str(r.get("Chromosome", "")).replace("chr", "").strip()
                pos = str(int(float(r.get("Start_Position", 0))))
                key = f"chr{ch}:{pos}"
                locus_meta[key] = (r.get("Target_Gene", "Unknown"), r.get("Element_Class", "Regulatory"))

    canonical_drivers = {
        "TERT", "GATA3", "FOXA1", "PLEKHS1", "DPH3", "OXNAD1", "RPL13A", "VRK3",
        "SDHD", "NEAT1", "MALAT1", "BCL2", "IGH", "PAX5", "BCL6", "NFKBIZ",
        "ZFP36L2", "ESR1", "TBC1D12", "MRPS31", "PES1", "ALDOA", "ST3GAL3", "AP3D1", "RAPSN"
    }

    with open(ml_path, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    X: list[list[float]] = []
    y: list[int] = []
    meta_rows: list[dict] = []

    for r in rows:
        loc = r["Locus"].strip()
        gene_str, elem_str = locus_meta.get(loc, ("Unknown", "Regulatory"))
        genes = {g.strip().upper() for g in gene_str.replace(",", ";").split(";") if g.strip()}
        is_known = 1 if (genes & canonical_drivers) else 0

        if not genes or genes == {"UNKNOWN"}:
            is_known = 1 if (
                int(float(r["Mutated_Donors"])) >= 5
                and int(float(r["Has_TF_Impact"])) == 1
                and (int(float(r["Is_Promoter"])) == 1 or int(float(r["Is_Enhancer"])) == 1)
            ) else 0

        feats = [
            float(r["Mutated_Donors"]),
            math.log10(max(float(r["Fold_Enrichment"]), 0.0) + 1.0),
            math.log10(max(float(r["Absolute_Distance_to_TSS"]), 0.0) + 1.0),
            float(r["Is_Promoter"]),
            float(r["Is_Enhancer"]),
            float(r["Has_TF_Impact"]),
            float(r["Motif_Gain"]),
            float(r["Motif_Loss"]),
        ]
        X.append(feats)
        y.append(is_known)
        meta_rows.append({
            "Locus": loc,
            "Cancer_Type": r["Cancer_Type"],
            "Target_Gene": gene_str,
            "Element_Class": elem_str,
            "Mutated_Donors": int(float(r["Mutated_Donors"])),
            "Fold_Enrichment": float(r["Fold_Enrichment"]),
            "Has_TF_Impact": int(float(r["Has_TF_Impact"])),
            "Benchmark_Label": is_known,
        })

    n = len(y)
    print(f"Benchmark Cohort: {n} noncoding hotspot loci (Positive={sum(y)}, Negative={n - sum(y)})")

    oof_gbdt = []
    oof_logreg = []
    for i in range(n):
        X_tr = [X[k] for k in range(n) if k != i]
        y_tr = [y[k] for k in range(n) if k != i]
        oof_gbdt.append(fit_predict_gbdt_stumps(X_tr, y_tr, X[i]))
        oof_logreg.append(fit_predict_logistic_l2(X_tr, y_tr, X[i]))

    base_fe = [float(r["Fold_Enrichment"]) for r in meta_rows]
    base_donors = [float(r["Mutated_Donors"]) for r in meta_rows]

    models = {
        "NC-DriverEngine GBDT (LOOCV)": (oof_gbdt, "#DD5143"),
        "L2 Logistic Regression (LOOCV)": (oof_logreg, "#4C72B0"),
        "Baseline: Fold Enrichment Only": (base_fe, "#2CA02C"),
        "Baseline: Donor Count Only": (base_donors, "#7F8C8D"),
    }

    curves_for_plot = {}
    metric_rows = []
    print("\nOut-of-Fold Cross-Validation Performance Summary:")
    print(f"{'Model / Ranking Criterion':34s} | {'ROC-AUC':8s} | {'PR-AUC (AP)':11s}")
    print("-" * 60)
    for m_name, (scores, col) in models.items():
        roc_a, pr_a, roc_pts, pr_pts = compute_roc_pr_curves(y, scores)
        curves_for_plot[m_name] = (roc_a, pr_a, roc_pts, pr_pts, col)
        metric_rows.append({
            "Model": m_name,
            "Evaluation": "Leave-One-Out CV (N=25)",
            "Positive_Loci": sum(y),
            "Negative_Loci": n - sum(y),
            "ROC_AUC": round(roc_a, 4),
            "PR_AUC": round(pr_a, 4),
        })
        print(f"{m_name:34s} | {roc_a:8.4f} | {pr_a:11.4f}")

    out_metrics = base / "results" / "tables" / "36_ml_loocv_benchmark_metrics.csv"
    with open(out_metrics, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(metric_rows[0].keys()))
        writer.writeheader()
        writer.writerows(metric_rows)

    for mrow, p_gb, p_lr in zip(meta_rows, oof_gbdt, oof_logreg):
        mrow["OOF_Prob_GBDT"] = round(p_gb, 4)
        mrow["OOF_Prob_LogReg"] = round(p_lr, 4)

    meta_rows.sort(key=lambda x: -x["OOF_Prob_GBDT"])
    out_preds = base / "results" / "tables" / "37_ml_loocv_locus_predictions.csv"
    with open(out_preds, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(meta_rows[0].keys()))
        writer.writeheader()
        writer.writerows(meta_rows)

    fig_dir = base / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    svg_out = fig_dir / "fig18_ml_pr_roc_benchmark.svg"
    _write_svg_pr_roc(svg_out, curves_for_plot)
    print(f"\nSaved LOOCV benchmark metrics  -> {out_metrics}")
    print(f"Saved per-locus OOF predictions -> {out_preds}")
    print(f"Saved PR/ROC vector SVG figure  -> {svg_out}")


if __name__ == "__main__":
    BASE_DIR = Path(r"D:\noncoding mutation finding surver")
    run_breast_nb_calibration(BASE_DIR)
    run_ml_loocv_benchmark(BASE_DIR)