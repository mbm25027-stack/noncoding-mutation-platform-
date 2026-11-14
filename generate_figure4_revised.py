import os
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.abspath(os.getcwd())
OUT_PNG = os.path.join(ROOT, "Figure4_ClinVar_Benchmark_Revised.png")
OUT_PDF = os.path.join(ROOT, "Figure4_ClinVar_Benchmark_Revised.pdf")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "axes.labelweight": "bold",
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 7.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def find_file(filename):
    direct = os.path.join(ROOT, filename)
    if os.path.isfile(direct):
        return direct
    for dp, dns, fns in os.walk(ROOT):
        dns[:] = [d for d in dns if d not in {"venv", ".venv", ".git", "__pycache__"}]
        if filename in fns:
            return os.path.join(dp, filename)
    return None


def count_vcf_pools(vcf_path):
    if not vcf_path or not os.path.isfile(vcf_path):
        return None, None, None
    chr_dev = {str(i) for i in range(1, 12)} | {"X"}
    chr_blind = {str(i) for i in range(12, 23)}
    n_dev = 0
    n_blind = 0
    tot = 0
    with open(vcf_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.split("\t", 1)[0].strip()
            cn = c[3:] if c.lower().startswith("chr") else c
            if cn in chr_dev:
                n_dev += 1
            elif cn in chr_blind:
                n_blind += 1
            tot += 1
    return tot, n_dev, n_blind


def check_frozen_splice_length(py_path):
    if not py_path or not os.path.isfile(py_path):
        return "Unknown (nc_mutscan_research_suite.py not found)"
    with open(py_path, "r", encoding="utf-8", errors="ignore") as f:
        txt = f.read()
    m = re.search(r"SPLICE_3SS.*?\[(.*?)\]", txt, re.DOTALL)
    if m:
        rows = re.findall(r"\{", m.group(1))
        if rows:
            return f"{len(rows)}-mer (detected in SPLICE_3SS)"
    for line in txt.splitlines():
        if "acceptor" in line.lower() and any(k in line for k in ("-20", "23", "-10", "13")):
            return line.strip()[:90]
    return "Check SPLICE_3SS definition manually"


def compute_metrics(tp, fp, tn, fn):
    sens = tp / (tp + fn) * 100.0
    spec = tn / (tn + fp) * 100.0
    ppv = tp / (tp + fp) * 100.0
    f1 = (2.0 * tp) / (2.0 * tp + fp + fn) * 100.0
    acc = (tp + tn) / (tp + fp + tn + fn) * 100.0
    fpr = fp / (tn + fp) * 100.0
    return {
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "sens": sens, "spec": spec, "ppv": ppv, "f1": f1, "acc": acc, "fpr": fpr,
    }


def main():
    clinvar_path = find_file("clinvar.vcf")
    scanner_path = find_file("nc_mutscan_research_suite.py")

    tot_vcf, dev_pool_n, blind_pool_n = count_vcf_pools(clinvar_path)
    if dev_pool_n is None:
        dev_pool_n = 955
    if blind_pool_n is None:
        blind_pool_n = 150

    splice_info = check_frozen_splice_length(scanner_path)
    print(f"[INFO] clinvar.vcf path      : {clinvar_path}")
    print(f"[INFO] chr1-11,X dev pool N  : {dev_pool_n}")
    print(f"[INFO] chr12-22 held-out X   : {blind_pool_n}  --> Replace [X] in Section 3.3 with {blind_pool_n}")
    print(f"[INFO] Frozen splice acceptor : {splice_info}")

    m_single_dev = compute_metrics(tp=50, fp=27, tn=48, fn=25)
    m_diff_dev = compute_metrics(tp=30, fp=8, tn=67, fn=45)
    m_diff_blind = compute_metrics(tp=37, fp=1, tn=74, fn=38)

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.4), gridspec_kw={"width_ratios": [1.05, 1.30, 0.95]})
    fig.suptitle(
        "Chromosome-Stratified ClinVar Benchmark of Frozen NC-MutScan v3.3",
        fontsize=12,
        fontweight="bold",
        y=1.02,
    )

    axa = axes[0]
    cats = ["chr1-11,X\nPool", "chr1-11,X\nDev Sample", "chr12-22\nHeld-Out", "chr12-22\nBlind Sample"]
    vals = [dev_pool_n, 150, blind_pool_n, 150]
    cols = ["#2b5c8f", "#5a9bd4", "#1b6a5a", "#4db09b"]
    bars_a = axa.bar(cats, vals, color=cols, width=0.55, edgecolor="#222222", linewidth=0.8)
    axa.set_yscale("log")
    axa.set_yticks([10, 100, 1000, 10000])
    axa.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    axa.set_ylim(10, max(vals) * 5.0)
    axa.set_ylabel("Variant Count (log scale)")
    axa.set_title("A. Chromosome Stratification", loc="left", pad=10)
    axa.grid(axis="y", linestyle="--", alpha=0.35)
    axa.set_axisbelow(True)
    axa.tick_params(axis="x", labelsize=7.8)

    annots_a = [
        f"N={dev_pool_n}\n(471 P, 484 B)",
        "N=150\n(75 P, 75 B)",
        f"N={blind_pool_n}\n(Held-Out Pool)",
        "N=150\n(75 P, 75 B)",
    ]
    for bar, txt in zip(bars_a, annots_a):
        h = bar.get_height()
        axa.text(
            bar.get_x() + bar.get_width() / 2.0,
            h * 1.18,
            txt,
            ha="center",
            va="bottom",
            fontsize=7.5,
            fontweight="bold",
            color="#222222",
        )

    axb = axes[1]
    metrics_labels = ["Sensitivity\n(TPR)", "Specificity\n(TNR)", "Precision\n(PPV)", "F1 Score\n(x100)"]
    x = np.arange(len(metrics_labels))
    w = 0.26
    s1 = [m_single_dev["sens"], m_single_dev["spec"], m_single_dev["ppv"], m_single_dev["f1"]]
    s2 = [m_diff_dev["sens"], m_diff_dev["spec"], m_diff_dev["ppv"], m_diff_dev["f1"]]
    s3 = [m_diff_blind["sens"], m_diff_blind["spec"], m_diff_blind["ppv"], m_diff_blind["f1"]]

    b1 = axb.bar(x - w, s1, width=w, label="Single-Seq PWM (Dev, N=150)", color="#a0aab4", edgecolor="#222222", linewidth=0.7)
    b2 = axb.bar(x, s2, width=w, label="Diff v3.3 (Dev Cal., N=150)", color="#3b75af", edgecolor="#222222", linewidth=0.7)
    b3 = axb.bar(x + w, s3, width=w, label="Diff v3.3 (Locked Blind, N=150)", color="#1b7a56", edgecolor="#222222", linewidth=0.7)

    axb.set_ylim(0, 122)
    axb.set_yticks([0, 20, 40, 60, 80, 100])
    axb.set_xticks(x)
    axb.set_xticklabels(metrics_labels)
    axb.set_ylabel("Percentage / Score (%)")
    axb.set_title("B. Dev vs. Locked Blind Performance", loc="left", pad=10)
    axb.grid(axis="y", linestyle="--", alpha=0.35)
    axb.set_axisbelow(True)
    axb.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92)

    for bars in (b1, b2, b3):
        for b in bars:
            h = b.get_height()
            axb.text(
                b.get_x() + b.get_width() / 2.0,
                h + 1.5,
                f"{h:.1f}",
                ha="center",
                va="bottom",
                fontsize=7.0,
                fontweight="bold",
                color="#222222",
            )

    axc = axes[2]
    cm = np.array([[37, 38], [1, 74]], dtype=float)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("cm_blues", ["#f4f7fa", "#1b5e91"])
    axc.imshow(cm, cmap=cmap, vmin=0, vmax=80)
    axc.spines["top"].set_visible(True)
    axc.spines["right"].set_visible(True)
    axc.set_xticks([0, 1])
    axc.set_yticks([0, 1])
    axc.set_xticklabels(["Predicted\nDriver (+)", "Predicted\nNeutral (-)"])
    axc.set_yticklabels(["ClinVar\nPathogenic", "ClinVar\nBenign"])
    axc.set_title("C. Locked Blind Matrix (N=150)", loc="left", pad=10)

    cell_labels = [
        [("TP = 37", "49.33% Sens"), ("FN = 38", "50.67% Miss")],
        [("FP = 1", "1.33% FPR"), ("TN = 74", "98.67% Spec")],
    ]
    for r in range(2):
        for c in range(2):
            val = cm[r, c]
            txt_col = "white" if val > 45 else "#111111"
            main_t, sub_t = cell_labels[r][c]
            axc.text(c, r - 0.08, main_t, ha="center", va="center", fontsize=10.0, fontweight="bold", color=txt_col)
            axc.text(c, r + 0.18, sub_t, ha="center", va="center", fontsize=8.0, color=txt_col)

    axc.set_xlabel("Acc: 74.00% | FPR: 1.33% | PPV: 97.37%", fontsize=8.2, fontweight="bold", labelpad=8)

    plt.tight_layout()
    fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    print(f"[SUCCESS] Saved Figure 4 PNG: {OUT_PNG}")
    print(f"[SUCCESS] Saved Figure 4 PDF: {OUT_PDF}")


if __name__ == "__main__":
    main()