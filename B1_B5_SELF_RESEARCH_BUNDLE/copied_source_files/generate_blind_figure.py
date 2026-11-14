import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.9

fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0), dpi=300, gridspec_kw={"width_ratios": [1.05, 1.45, 1.05]})
plt.subplots_adjust(wspace=0.34, top=0.82, bottom=0.18, left=0.06, right=0.98)

# ------------------------------------------------------------------------------
# PANEL A: Zero-Leakage Dataset Construction
# ------------------------------------------------------------------------------
ax = axes[0]
categories = ["Dev Pool\n(chr1–11)", "Blind Pool\n(chr12–22)", "Dev Sample\n(chr1–11)", "Blind Holdout\n(chr12–22)"]
path_counts = [6036, 3472, 75, 75]
benign_counts = [314319, 193490, 75, 75]

y_pos = np.arange(len(categories))
ax.barh(y_pos - 0.18, path_counts, height=0.34, color="#C0392B", label="Pathogenic (Non-Coding)")
ax.barh(y_pos + 0.18, benign_counts, height=0.34, color="#2E86C1", label="Benign (Non-Coding)")
ax.set_xscale("log")
ax.set_xlim(10, 2_000_000)
ax.set_yticks(y_pos)
ax.set_yticklabels(categories, fontsize=9.5, fontweight="bold")
ax.invert_yaxis()
ax.set_xlabel("Variant Count (Log Scale)", fontsize=9.5, fontweight="bold")
ax.set_title("A. Isoform-Purged Chromosome Splits", fontsize=11, fontweight="bold", loc="left", pad=10)
ax.grid(axis="x", linestyle="--", alpha=0.35)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for idx, (p, b) in enumerate(zip(path_counts, benign_counts)):
    ax.text(p * 1.25, idx - 0.18, f"{p:,}", va="center", fontsize=8, color="#922B21", fontweight="bold")
    ax.text(b * 1.25, idx + 0.18, f"{b:,}", va="center", fontsize=8, color="#1B4F72", fontweight="bold")

ax.legend(loc="lower right", frameon=True, fontsize=8)

# ------------------------------------------------------------------------------
# PANEL B: Version Ablation & Blind Generalization
# ------------------------------------------------------------------------------
ax = axes[1]
metrics = ["Sensitivity\n(TPR)", "Specificity\n(TNR)", "Precision\n(PPV)", "F1-Score\n(×100)"]
v31_dev = [72.00, 21.30, 47.79, 57.45]
v33_dev = [66.67, 89.33, 86.21, 75.19]
v33_blind = [49.33, 98.67, 97.37, 65.49]

x = np.arange(len(metrics))
w = 0.26

b1 = ax.bar(x - w, v31_dev, width=w, color="#95A5A6", label="v3.1 Single-Seq (Dev)")
b2 = ax.bar(x, v33_dev, width=w, color="#2980B9", label="v3.3 ΔPWM (Dev chr1–11)")
b3 = ax.bar(x + w, v33_blind, width=w, color="#1E8449", label="v3.3 ΔPWM (Blind chr12–22)")

ax.set_ylim(0, 114)
ax.set_xticks(x)
ax.set_xticklabels(metrics, fontsize=9.5, fontweight="bold")
ax.set_ylabel("Percentage / Score (%)", fontsize=9.5, fontweight="bold")
ax.set_title("B. Differential ΔPWM Specificity Gain", fontsize=11, fontweight="bold", loc="left", pad=10)
ax.grid(axis="y", linestyle="--", alpha=0.35)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

for bars in (b1, b2, b3):
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 1.5, f"{h:.1f}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")

ax.legend(loc="upper right", frameon=True, fontsize=8)

# ------------------------------------------------------------------------------
# PANEL C: Locked Blind Holdout Confusion Matrix (chr12-22, n=150)
# ------------------------------------------------------------------------------
ax = axes[2]
cm = np.array([[37, 38], [1, 74]])
cm_pct = np.array([[37 / 75 * 100, 38 / 75 * 100], [1 / 75 * 100, 74 / 75 * 100]])
cell_labels = [["TP", "FN"], ["FP", "TN"]]

im = ax.imshow(cm, cmap="Blues", vmin=0, vmax=80)
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(["Predicted\nDriver", "Predicted\nBenign"], fontsize=9.5, fontweight="bold")
ax.set_yticklabels(["ClinVar\nPathogenic", "ClinVar\nBenign"], fontsize=9.5, fontweight="bold")
ax.set_title("C. Locked Blind Matrix (chr12–22)", fontsize=11, fontweight="bold", loc="left", pad=10)

for i in range(2):
    for j in range(2):
        txt_color = "white" if cm[i, j] > 45 else "#111111"
        ax.text(
            j, i,
            f"{cell_labels[i][j]} = {cm[i, j]}\n({cm_pct[i, j]:.1f}%)",
            ha="center", va="center", color=txt_color, fontsize=10, fontweight="bold"
        )

ax.text(
    0.5, -0.24,
    "Accuracy: 74.00%  |  FPR: 1.33%  |  PPV: 97.37%",
    transform=ax.transAxes, ha="center", fontsize=8.5, fontweight="bold", color="#1E8449"
)

fig.suptitle("NC-MutScan v3.3 Zero-Leakage Chromosome-Stratified Benchmark Evaluation", fontsize=13, fontweight="bold", y=0.96)
plt.savefig("Figure_19_Blind_Benchmark_v33.png", dpi=300)
plt.savefig("Figure_19_Blind_Benchmark_v33.pdf")
print("Saved Figure_19_Blind_Benchmark_v33.png and Figure_19_Blind_Benchmark_v33.pdf")
