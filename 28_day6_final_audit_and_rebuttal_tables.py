import csv
import math
import random
from pathlib import Path


def compute_roc_auc(y_true: list[int], y_score: list[float]) -> float:
    pos = [s for y, s in zip(y_true, y_score) if y == 1]
    neg = [s for y, s in zip(y_true, y_score) if y == 0]
    if not pos or not neg:
        return 0.5
    wins = sum(1.0 if sp > sn else (0.5 if sp == sn else 0.0) for sp in pos for sn in neg)
    return wins / (len(pos) * len(neg))


def compute_pr_auc(y_true: list[int], y_score: list[float]) -> float:
    n_pos = sum(y_true)
    if n_pos == 0:
        return 0.0
    paired = sorted(zip(y_score, y_true), key=lambda x: -x[0])
    tp = 0
    ap = 0.0
    for rank, (_, y) in enumerate(paired, start=1):
        if y == 1:
            tp += 1
            ap += tp / float(rank)
    return ap / float(n_pos)


def sigmoid(z: float) -> float:
    return 1.0 / (1.0 + math.exp(-max(min(z, 60.0), -60.0)))


def fit_predict_l2_logreg(X_tr: list[list[float]], y_tr: list[int], x_te: list[float], l2: float = 0.4) -> float:
    n, d = len(X_tr), len(X_tr[0])
    w, b, lr = [0.0] * d, 0.0, 0.15
    for _ in range(350):
        gw, gb = [0.0] * d, 0.0
        for xi, yi in zip(X_tr, y_tr):
            err = sigmoid(b + sum(wj * xij for wj, xij in zip(w, xi))) - float(yi)
            gb += err
            for j in range(d):
                gw[j] += err * xi[j]
        for j in range(d):
            w[j] -= lr * (gw[j] / n + l2 * w[j])
        b -= lr * (gb / n)
    return sigmoid(b + sum(wj * xj for wj, xj in zip(w, x_te)))


def fit_predict_gbdt(
    X_tr: list[list[float]], y_tr: list[int], x_te: list[float],
    n_estimators: int = 20, lr: float = 0.18, reg_lambda: float = 1.2
) -> float:
    n, d = len(X_tr), len(X_tr[0])
    p0 = min(max(sum(y_tr) / float(n), 1e-4), 1.0 - 1e-4)
    f0 = math.log(p0 / (1.0 - p0))
    F = [f0] * n
    stumps = []
    for _ in range(n_estimators):
        probs = [sigmoid(fi) for fi in F]
        g = [pi - float(yi) for pi, yi in zip(probs, y_tr)]
        h = [max(pi * (1.0 - pi), 1e-6) for pi in probs]
        best_gain = -1e18
        best_stump = (0, 0.0, 0.0, 0.0)
        G_tot, H_tot = sum(g), sum(h)
        for j in range(d):
            vals = sorted(set(row[j] for row in X_tr))
            if len(vals) < 2:
                continue
            for k in range(len(vals) - 1):
                thr = 0.5 * (vals[k] + vals[k + 1])
                GL = sum(g[i] for i in range(n) if X_tr[i][j] <= thr)
                HL = sum(h[i] for i in range(n) if X_tr[i][j] <= thr)
                GR, HR = G_tot - GL, H_tot - HL
                gain = 0.5 * (
                    (GL * GL) / (HL + reg_lambda)
                    + (GR * GR) / (HR + reg_lambda)
                    - (G_tot * G_tot) / (H_tot + reg_lambda)
                )
                if gain > best_gain:
                    best_gain = gain
                    best_stump = (j, thr, -GL / (HL + reg_lambda), -GR / (HR + reg_lambda))
        jb, tb, wL, wR = best_stump
        stumps.append(best_stump)
        for i in range(n):
            F[i] += lr * (wL if X_tr[i][jb] <= tb else wR)
    f_te = f0 + sum(lr * (wL if x_te[jb] <= tb else wR) for jb, tb, wL, wR in stumps)
    return sigmoid(f_te)


def run_final_loocv_ci(base: Path):
    tab = base / "results" / "tables"

    # 1. Synchronize Table 42 directly from 31_covariate_adjusted_eqtl_results.csv
    eqtl_31 = list(csv.DictReader(open(tab / "31_covariate_adjusted_eqtl_results.csv", "r", encoding="utf-8")))
    order = {"BCL6": 0, "BCL2": 1, "PAX5": 2, "WDR74": 3}
    eqtl_31.sort(key=lambda r: order.get(r["Gene"], 99))
    out_42 = []
    for r in eqtl_31:
        out_42.append({
            "Gene": r["Gene"],
            "Dataset": r["Dataset"],
            "Locus_hg19": r["Locus_GRCh38"],
            "CompleteCase_WT_vs_MUT": f"{r['WT_Donors']} / {r['Mut_Donors']}",
            "Mean_WT_vs_MUT": f"{float(r['Mean_WT_Expr']):.3f} vs. {float(r['Mean_Mut_Expr']):.3f}",
            "Unadjusted_P": f"{float(r['Unadj_P']):.2e}" if float(r["Unadj_P"]) < 0.01 else f"{float(r['Unadj_P']):.3f}",
            "CovAdj_Beta_Mut": f"{float(r['CovAdj_Beta_Mut']):+.3f}",
            "CovAdj_P": f"{float(r['CovAdj_P']):.2e}" if float(r["CovAdj_P"]) < 0.01 else f"{float(r['CovAdj_P']):.3f}",
            "Q_eQTL_BH": f"{float(r['Q_eQTL_BH']):.2e}" if float(r["Q_eQTL_BH"]) < 0.01 else f"{float(r['Q_eQTL_BH']):.3f}",
        })
    with open(tab / "42_eqtl_unadjusted_and_complete_case.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_42[0].keys()))
        w.writeheader()
        w.writerows(out_42)

    # 2. Exact LOOCV + 2,000 Stratified Bootstrap 95% CIs + 2,000 Label Permutations
    mat_rows = list(csv.DictReader(open(tab / "38_repaired_ml_hotspot_matrix.csv", "r", encoding="utf-8")))
    n = len(mat_rows)
    loci, cancers, y_drv, y_arch, X_drv, X_arch = [], [], [], [], [], []
    for r in mat_rows:
        ca = r["Cancer_Type"].strip()
        is_mel = 1 if "mel" in ca.lower() else 0
        drv = int(float(r["NB_Significant_Driver"]))
        log_fe = math.log10(max(float(r["Fold_Enrichment_Repaired"]), 0.0) + 1.0)
        log_tss = math.log10(abs(float(r["Absolute_Distance_to_TSS"])) + 1.0)
        is_prom = float(r["Is_Promoter"])
        is_enh = float(r["Is_Enhancer"])
        has_tf = float(r["Has_TF_Impact"])
        m_gain = float(r["Motif_Gain"])
        m_loss = float(r["Motif_Loss"])

        loci.append(r["Locus"].strip())
        cancers.append(ca)
        y_drv.append(drv)
        y_arch.append(is_mel)
        X_drv.append([log_fe, log_tss, is_prom, is_enh, has_tf, m_gain, m_loss])
        X_arch.append([log_tss, is_prom, is_enh, has_tf, m_gain, m_loss])

    def eval_loocv(X_raw: list[list[float]], y: list[int], baselines: dict[str, list[float]]):
        d = len(X_raw[0])
        oof_g, oof_l = [], []
        for i in range(n):
            Xtr = [X_raw[k] for k in range(n) if k != i]
            ytr = [y[k] for k in range(n) if k != i]
            xte = X_raw[i]
            mu = [sum(row[j] for row in Xtr) / len(Xtr) for j in range(d)]
            sd = [max(math.sqrt(sum((row[j] - mu[j]) ** 2 for row in Xtr) / len(Xtr)), 1e-8) for j in range(d)]
            Xtr_s = [[(row[j] - mu[j]) / sd[j] for j in range(d)] for row in Xtr]
            xte_s = [(xte[j] - mu[j]) / sd[j] for j in range(d)]
            oof_g.append(fit_predict_gbdt(Xtr_s, ytr, xte_s))
            oof_l.append(fit_predict_l2_logreg(Xtr_s, ytr, xte_s))
        res = {"NC-DriverEngine GBDT (LOOCV)": oof_g, "L2 Logistic Regression (LOOCV)": oof_l}
        res.update(baselines)
        return res

    s1 = eval_loocv(
        X_drv, y_drv,
        {"Baseline: TSS Proximity Only": [-r[1] for r in X_drv], "Baseline: Fold Enrichment Only": [r[0] for r in X_drv]}
    )
    s2 = eval_loocv(
        X_arch, y_arch,
        {"Baseline: TSS Proximity Only": [-r[0] for r in X_arch], "Baseline: TF Motif Impact Only": [r[3] for r in X_arch]}
    )

    pred_rows = []
    for tname, y_vec, sdict in [("NB_Driver_Prioritization", y_drv, s1), ("Tumor_Regulatory_Architecture", y_arch, s2)]:
        for i in range(n):
            pr = {"Benchmark_Task": tname, "Locus": loci[i], "Cancer": cancers[i], "True_Label": y_vec[i]}
            for mname, svec in sdict.items():
                pr[mname] = round(svec[i], 6)
            pred_rows.append(pr)
    with open(tab / "37_loocv_predictions.csv", "w", newline="", encoding="utf-8") as f:
        fields = list(dict.fromkeys(k for r in pred_rows for k in r.keys()))
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(pred_rows)

    rng = random.Random(20261113)
    out_43 = []
    for tname, y_true, sdict in [("NB_Driver_Prioritization", y_drv, s1), ("Tumor_Regulatory_Architecture", y_arch, s2)]:
        p_idx = [i for i, y in enumerate(y_true) if y == 1]
        n_idx = [i for i, y in enumerate(y_true) if y == 0]
        for mname, sc in sdict.items():
            roc = compute_roc_auc(y_true, sc)
            pr = compute_pr_auc(y_true, sc)
            b_roc, b_pr = [], []
            for _ in range(2000):
                bi = [rng.choice(p_idx) for _ in p_idx] + [rng.choice(n_idx) for _ in n_idx]
                by = [y_true[i] for i in bi]
                bs = [sc[i] for i in bi]
                b_roc.append(compute_roc_auc(by, bs))
                b_pr.append(compute_pr_auc(by, bs))
            b_roc.sort()
            b_pr.sort()
            pr_p = roc_p = 0
            yp = list(y_true)
            for _ in range(2000):
                rng.shuffle(yp)
                if compute_roc_auc(yp, sc) >= roc:
                    roc_p += 1
                if compute_pr_auc(yp, sc) >= pr:
                    pr_p += 1
            out_43.append({
                "Benchmark_Task": tname,
                "Model": mname,
                "LOOCV_ROC_AUC": round(roc, 4),
                "ROC_95CI": f"[{b_roc[50]:.4f}, {b_roc[1949]:.4f}]",
                "Perm_P_ROC": f"{(roc_p + 1) / 2001.0:.4f}",
                "LOOCV_PR_AUC": round(pr, 4),
                "PR_95CI": f"[{b_pr[50]:.4f}, {b_pr[1949]:.4f}]",
                "Perm_P_PR": f"{(pr_p + 1) / 2001.0:.4f}",
            })

    with open(tab / "43_loocv_bootstrap_and_permutation.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_43[0].keys()))
        w.writeheader()
        w.writerows(out_43)

    print("=== FINAL TABLE 8 (42_eqtl_unadjusted_and_complete_case.csv) ===")
    for r in out_42:
        print(f"  {r['Gene']:6s} | WT/MUT={r['CompleteCase_WT_vs_MUT']:8s} | Unadj_P={r['Unadjusted_P']:8s} | Beta={r['CovAdj_Beta_Mut']} | Q_eQTL={r['Q_eQTL_BH']}")

    print("\n=== FINAL LOOCV + 95% CI + PERMUTATION P (43_loocv_bootstrap_and_permutation.csv) ===")
    for s in out_43:
        print(
            f"  [{s['Benchmark_Task'][:22]:22s}] {s['Model']:30s} | "
            f"PR-AUC={s['LOOCV_PR_AUC']:.4f} {s['PR_95CI']} (p_perm={s['Perm_P_PR']}) | "
            f"ROC-AUC={s['LOOCV_ROC_AUC']:.4f} {s['ROC_95CI']}"
        )


if __name__ == "__main__":
    run_final_loocv_ci(Path(r"D:\noncoding mutation finding surver"))