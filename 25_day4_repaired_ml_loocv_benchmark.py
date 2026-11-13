import csv
import math
from pathlib import Path


# ==============================================================================
# 1. Pure-Python Metrics & Leakage-Free LOOCV Models (Zero External DLLs)
# ==============================================================================

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
    l2: float = 0.4,
    lr: float = 0.15,
    epochs: int = 450,
) -> float:
    """In-fold standardized L2 Logistic Regression (strictly zero test leakage)."""
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
            z = max(min(b + sum(w[j] * Xs[i][j] for j in range(d)), 30.0), -30.0)
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
    n_estimators: int = 20,
    lr: float = 0.18,
    lam: float = 1.2,
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
# 2. Repair Missing Melanoma Fold_Enrichment & Join Calibrated NB Q-Values
# ==============================================================================

def load_and_repair_matrix(base: Path) -> list[dict]:
    ml_path = base / "results" / "tables" / "ml_hotspot_training_matrix.csv"
    mel_nb_path = base / "results" / "tables" / "33_melanoma_recalibrated_nb_clusters.csv"
    mel_cand_path = base / "results" / "tables" / "26_final_noncoding_driver_candidates.csv"
    brc_nb_path = base / "results" / "tables" / "35_breast_recalibrated_candidates.csv"
    out_repaired = base / "results" / "tables" / "38_repaired_ml_hotspot_matrix.csv"

    mel_records = []
    if mel_nb_path.exists():
        with open(mel_nb_path, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                ch = str(r.get("Chromosome", "")).replace("chr", "").strip()
                pos = int(float(r.get("Peak_Center", 0)))
                fe = float(r.get("Local_Obs_Over_Exp", 0.0))
                gene = r.get("Target_Gene", "Unknown")
                q_nb = float(r.get("Q_NB_ZeroInc", 1.0))
                mel_records.append((ch, pos, fe, gene, q_nb))
    elif mel_cand_path.exists():
        with open(mel_cand_path, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                ch = str(r.get("Chromosome", "")).replace("chr", "").strip()
                pos = int(float(r.get("Center_Position", 0)))
                fe = float(r.get("Local_Obs_Over_Exp") or r.get("Obs_Over_Exp") or 0.0)
                gene = r.get("Target_Gene", "Unknown")
                q_nb = float(r.get("q_value_local") or r.get("q_value") or 1.0)
                mel_records.append((ch, pos, fe, gene, q_nb))

    brc_records = []
    if brc_nb_path.exists():
        with open(brc_nb_path, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                ch = str(r.get("Chromosome", "")).replace("chr", "").strip()
                pos = int(float(r.get("Start_Position", 0)))
                fe = float(r.get("Fold_Enrichment", 0.0))
                gene = r.get("Target_Gene", "Unknown")
                q_nb = float(r.get("Q_NB_ZeroInc", 1.0))
                brc_records.append((ch, pos, fe, gene, q_nb))

    with open(ml_path, "r", encoding="utf-8") as f:
        raw_rows = list(csv.DictReader(f))

    repaired_rows = []
    for r in raw_rows:
        loc = r["Locus"].strip()
        ch_part, pos_part = loc.split(":")
        ch = ch_part.replace("chr", "").strip()
        pos = int(float(pos_part))
        ctype = r["Cancer_Type"].strip()
        fe = float(r["Fold_Enrichment"])
        gene = "Unknown"
        q_nb = 1.0

        pool = mel_records if ctype == "Melanoma" else brc_records
        best_hit = None
        best_dist = 10**9
        for r_ch, r_pos, r_fe, r_gene, r_q in pool:
            if r_ch == ch and abs(r_pos - pos) < best_dist:
                best_dist = abs(r_pos - pos)
                best_hit = (r_fe, r_gene, r_q)

        if best_hit is not None and best_dist <= 250:
            if fe <= 0.0 and best_hit[0] > 0.0:
                fe = best_hit[0]
            gene = best_hit[1]
            q_nb = best_hit[2]

        if fe <= 0.0:
            fe = round(float(r["Mutated_Donors"]) / 0.12909, 2)

        is_nb_driver = 1 if (q_nb < 0.05) else 0

        repaired_rows.append({
            "Locus": loc,
            "Cancer_Type": ctype,
            "Target_Gene": gene,
            "Mutated_Donors": int(float(r["Mutated_Donors"])),
            "Fold_Enrichment_Repaired": round(fe, 2),
            "Distance_to_TSS": int(float(r["Distance_to_TSS"])),
            "Absolute_Distance_to_TSS": int(float(r["Absolute_Distance_to_TSS"])),
            "Is_Promoter": int(float(r["Is_Promoter"])),
            "Is_Enhancer": int(float(r["Is_Enhancer"])),
            "Has_TF_Impact": int(float(r["Has_TF_Impact"])),
            "Motif_Gain": int(float(r["Motif_Gain"])),
            "Motif_Loss": int(float(r["Motif_Loss"])),
            "Q_NB_ZeroInc": f"{q_nb:.3e}",
            "NB_Significant_Driver": is_nb_driver,
        })

    with open(out_repaired, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(repaired_rows[0].keys()))
        writer.writeheader()
        writer.writerows(repaired_rows)
    print(f"Saved repaired ML feature matrix (zero missing Fold_Enrichment) -> {out_repaired}")
    return repaired_rows


# ==============================================================================
# 3. Publication Vector SVG Generator (2-Panel ROC & PR Curves)
# ==============================================================================

def write_svg_pr_roc(
    svg_path: Path,
    title_main: str,
    curves: dict[str, tuple[float, float, list[tuple[float, float]], list[tuple[float, float]], str]],
):
    W, H = 1080, 480
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" style="background:#ffffff;font-family:Arial,Helvetica,sans-serif;">',
        f'<text x="{W/2}" y="30" text-anchor="middle" font-size="16" font-weight="bold" fill="#111111">{title_main}</text>',
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
            lx = ox + 105 if is_roc else ox + 15
            metric_txt = f"ROC-AUC = {roc_auc:.3f}" if is_roc else f"PR-AUC = {pr_auc:.3f}"
            lines.append(f'<line x1="{lx}" y1="{ly-4}" x2="{lx+22}" y2="{ly-4}" stroke="{col}" stroke-width="2.4"{dash}/>')
            lines.append(f'<text x="{lx+28}" y="{ly}" font-size="10.2" fill="#111111">{m_name} ({metric_txt})</text>')

    lines.append("</svg>")
    svg_path.write_text("\n".join(lines), encoding="utf-8")


# ==============================================================================
# 4. Run Dual LOOCV Benchmarks (Driver Prioritization & Regulatory Topology)
# ==============================================================================

def evaluate_loocv_suite(base: Path, rows: list[dict]):
    n = len(rows)
    fig_dir = base / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    out_summary = base / "results" / "tables" / "36_ml_loocv_benchmark_metrics.csv"
    out_preds = base / "results" / "tables" / "37_ml_loocv_locus_predictions.csv"
    all_metrics = []

    # --------------------------------------------------------------------------
    # TASK 1: Genome-Wide NB Driver Prioritization (Q_NB < 0.05 vs Sub-threshold)
    # --------------------------------------------------------------------------
    y_drv = [r["NB_Significant_Driver"] for r in rows]
    X_drv = [
        [
            math.log10(r["Fold_Enrichment_Repaired"] + 1.0),
            math.log10(r["Absolute_Distance_to_TSS"] + 1.0),
            float(r["Is_Promoter"]),
            float(r["Is_Enhancer"]),
            float(r["Has_TF_Impact"]),
            float(r["Motif_Gain"]),
            float(r["Motif_Loss"]),
        ]
        for r in rows
    ]

    oof_gbdt_drv = [
        fit_predict_gbdt_stumps([X_drv[k] for k in range(n) if k != i], [y_drv[k] for k in range(n) if k != i], X_drv[i])
        for i in range(n)
    ]
    oof_lr_drv = [
        fit_predict_logistic_l2([X_drv[k] for k in range(n) if k != i], [y_drv[k] for k in range(n) if k != i], X_drv[i])
        for i in range(n)
    ]
    base_fe = [r["Fold_Enrichment_Repaired"] for r in rows]
    base_tss = [-math.log10(r["Absolute_Distance_to_TSS"] + 1.0) for r in rows]

    models_drv = {
        "NC-DriverEngine GBDT (LOOCV)": (oof_gbdt_drv, "#DD5143"),
        "L2 Logistic Regression (LOOCV)": (oof_lr_drv, "#4C72B0"),
        "Baseline: TSS Proximity Only": (base_tss, "#2CA02C"),
        "Baseline: Fold Enrichment Only": (base_fe, "#7F8C8D"),
    }

    print(f"\n=== BENCHMARK 1: NB-Calibrated Driver Prioritization (Pos={sum(y_drv)}, Neg={n - sum(y_drv)}) ===")
    print(f"{'Model / Ranking Criterion':34s} | {'ROC-AUC':8s} | {'PR-AUC (AP)':11s}")
    print("-" * 60)
    curves_drv = {}
    for m_name, (scores, col) in models_drv.items():
        roc_a, pr_a, rpts, ppts = compute_roc_pr_curves(y_drv, scores)
        curves_drv[m_name] = (roc_a, pr_a, rpts, ppts, col)
        all_metrics.append({
            "Benchmark_Task": "NB_Driver_Prioritization (Q_NB < 0.05)",
            "Model": m_name,
            "Evaluation": "Strict Leave-One-Out CV (N=25)",
            "Positive_Class": sum(y_drv),
            "Negative_Class": n - sum(y_drv),
            "ROC_AUC": round(roc_a, 4),
            "PR_AUC": round(pr_a, 4),
        })
        print(f"{m_name:34s} | {roc_a:8.4f} | {pr_a:11.4f}")

    write_svg_pr_roc(
        fig_dir / "fig18_ml_pr_roc_benchmark.svg",
        f"Leakage-Free LOOCV Driver Prioritization Benchmark (N = 25 Loci; Pos = {sum(y_drv)}, Neg = {n - sum(y_drv)})",
        curves_drv,
    )

    # --------------------------------------------------------------------------
    # TASK 2: Tumor-Specific Regulatory Architecture (Melanoma=1 vs Breast=0)
    # --------------------------------------------------------------------------
    y_can = [1 if r["Cancer_Type"] == "Melanoma" else 0 for r in rows]
    X_can = [
        [
            math.log10(r["Absolute_Distance_to_TSS"] + 1.0),
            float(r["Is_Promoter"]),
            float(r["Is_Enhancer"]),
            float(r["Has_TF_Impact"]),
            float(r["Motif_Gain"]),
            float(r["Motif_Loss"]),
        ]
        for r in rows
    ]

    oof_gbdt_can = [
        fit_predict_gbdt_stumps([X_can[k] for k in range(n) if k != i], [y_can[k] for k in range(n) if k != i], X_can[i])
        for i in range(n)
    ]
    oof_lr_can = [
        fit_predict_logistic_l2([X_can[k] for k in range(n) if k != i], [y_can[k] for k in range(n) if k != i], X_can[i])
        for i in range(n)
    ]

    models_can = {
        "NC-DriverEngine GBDT (LOOCV)": (oof_gbdt_can, "#DD5143"),
        "L2 Logistic Regression (LOOCV)": (oof_lr_can, "#4C72B0"),
        "Baseline: TSS Proximity Only": (base_tss, "#2CA02C"),
        "Baseline: TF Motif Impact Only": ([float(r["Has_TF_Impact"]) for r in rows], "#7F8C8D"),
    }

    print(f"\n=== BENCHMARK 2: Scale-Invariant Regulatory Architecture (Melanoma={sum(y_can)} vs Breast={n - sum(y_can)}) ===")
    print(f"{'Model / Ranking Criterion':34s} | {'ROC-AUC':8s} | {'PR-AUC (AP)':11s}")
    print("-" * 60)
    curves_can = {}
    for m_name, (scores, col) in models_can.items():
        roc_a, pr_a, rpts, ppts = compute_roc_pr_curves(y_can, scores)
        curves_can[m_name] = (roc_a, pr_a, rpts, ppts, col)
        all_metrics.append({
            "Benchmark_Task": "Tumor_Regulatory_Architecture (Melanoma vs Breast)",
            "Model": m_name,
            "Evaluation": "Strict Leave-One-Out CV (N=25)",
            "Positive_Class": sum(y_can),
            "Negative_Class": n - sum(y_can),
            "ROC_AUC": round(roc_a, 4),
            "PR_AUC": round(pr_a, 4),
        })
        print(f"{m_name:34s} | {roc_a:8.4f} | {pr_a:11.4f}")

    write_svg_pr_roc(
        fig_dir / "fig19_regulatory_architecture_loocv.svg",
        "Scale-Invariant Regulatory Architecture Classification: Melanoma vs. Breast (LOOCV, N = 25)",
        curves_can,
    )

    with open(out_summary, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_metrics[0].keys()))
        writer.writeheader()
        writer.writerows(all_metrics)

    for r, p_d_gb, p_d_lr, p_c_gb, p_c_lr in zip(rows, oof_gbdt_drv, oof_lr_drv, oof_gbdt_can, oof_lr_can):
        r["OOF_Driver_Prob_GBDT"] = round(p_d_gb, 4)
        r["OOF_Driver_Prob_LogReg"] = round(p_d_lr, 4)
        r["OOF_Melanoma_Arch_Prob_GBDT"] = round(p_c_gb, 4)
        r["OOF_Melanoma_Arch_Prob_LogReg"] = round(p_c_lr, 4)

    with open(out_preds, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved combined LOOCV benchmark table  -> {out_summary}")
    print(f"Saved per-locus OOF predictions       -> {out_preds}")
    print(f"Saved Driver PR/ROC vector SVG        -> {fig_dir / 'fig18_ml_pr_roc_benchmark.svg'}")
    print(f"Saved Architecture PR/ROC vector SVG  -> {fig_dir / 'fig19_regulatory_architecture_loocv.svg'}")


if __name__ == "__main__":
    BASE_DIR = Path(r"D:\noncoding mutation finding surver")
    repaired = load_and_repair_matrix(BASE_DIR)
    evaluate_loocv_suite(BASE_DIR, repaired)