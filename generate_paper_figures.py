import matplotlib.pyplot as plt
import numpy as np

# Set publication-quality figure aesthetics
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 300,
    'font.family': 'sans-serif'
})

def generate_figure_1():
    """Figure 1: Multi-TF Competition Dynamics (Overlapping AP-1 / E2F1 Site)"""
    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    tfs = ['AP-1 (MA0099.2)\n[Secondary Gain]', 'E2F1 (MA0024.3)\n[Primary Driver Rank 1]']
    ref_scores = [3.0, 3.5]
    alt_scores = [7.2, 9.5]

    y = np.arange(len(tfs))
    height = 0.35

    ax.barh(y - height/2, ref_scores, height, label='Wild-Type (Ref) Score', color='#2b5c8f', edgecolor='black', linewidth=0.8)
    ax.barh(y + height/2, alt_scores, height, label='Mutant (Alt) Score', color='#c0392b', edgecolor='black', linewidth=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels(tfs)
    ax.set_xlabel('JASPAR PWM Log-Odds Score (bits)')
    ax.set_title('Figure 1: Multi-TF Competition Dynamics at Overlapping cis-Regulatory Locus', pad=15)
    ax.legend(loc='lower right', frameon=True)

    # Annotate delta shifts
    ax.text(7.4, y[0] + height/2, 'ΔPWM: +4.2 bits', va='center', color='#27ae60', fontweight='bold')
    ax.text(9.7, y[1] + height/2, 'ΔPWM: +6.0 bits', va='center', color='#27ae60', fontweight='bold')

    ax.set_xlim(0, 12)
    plt.tight_layout()
    plt.savefig('Fig1_Multi_TF_Competition.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('Fig1_Multi_TF_Competition.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()
    print("[✓] Generated Fig1_Multi_TF_Competition")

def generate_figure_2():
    """Figure 2: Dual Mechanisms for False-Positive Suppression (AID/SHM & GC-Noise)"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel A: AID/SHM Filter
    labels_a = ['Raw Matrix Score\n(Unfiltered)', 'Adjusted Score\n(AID RGYW Flag)']
    scores_a = [7.5, 3.2]
    bars_a = ax1.bar(labels_a, scores_a, color=['#e74c3c', '#f39c12'], width=0.45, edgecolor='black', linewidth=0.8)
    ax1.set_ylabel('Binding Affinity Score (bits)')
    ax1.set_title('A: AID/SHM Passenger Filter\n(Reverse-Strand RGYW Signature)', fontsize=12)
    ax1.set_ylim(0, 10)
    for bar in bars_a:
        yval = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2, yval + 0.3, f'{yval} bits', ha='center', va='bottom', fontweight='bold')

    # Panel B: High-GC Dirichlet Normalization
    labels_b = ['Legacy Scoring\n(False Positive)', 'NC-MutScan Engine\n(Dirichlet bg=0.25)']
    scores_b = [6.8, 0.0]
    bars_b = ax2.bar(labels_b, scores_b, color=['#7f8c8d', '#27ae60'], width=0.45, edgecolor='black', linewidth=0.8)
    ax2.set_ylabel('Binding Affinity Score (bits)')
    ax2.set_title('B: GC-Rich / G4 Noise Suppression\n(Background Normalization)', fontsize=12)
    ax2.set_ylim(0, 10)
    for bar in bars_b:
        yval = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2, yval + 0.3, f'{yval} bits', ha='center', va='bottom', fontweight='bold')

    plt.tight_layout()
    plt.savefig('Fig2_False_Positive_Suppression.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('Fig2_False_Positive_Suppression.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()
    print("[✓] Generated Fig2_False_Positive_Suppression")

def generate_figure_3():
    """Figure 3: Comparative Diagnostic Performance (NC-MutScan vs Legacy Tools)"""
    fig, ax = plt.subplots(figsize=(9, 5))

    tools = ['NC-MutScan v3.1', 'CADD v1.6', 'GWAVA', 'FunSeq2']
    sensitivity = [100.0, 64.2, 58.5, 61.0]
    specificity = [98.5, 71.0, 75.4, 68.2]

    x = np.arange(len(tools))
    width = 0.35

    rects1 = ax.bar(x - width/2, sensitivity, width, label='Sensitivity (%)', color='#2980b9', edgecolor='black', linewidth=0.8)
    rects2 = ax.bar(x + width/2, specificity, width, label='Specificity (%)', color='#e67e22', edgecolor='black', linewidth=0.8)

    ax.set_ylabel('Diagnostic Metric (%)')
    ax.set_title('Figure 3: Pan-Cancer Noncoding Hotspot Classification Performance', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(tools)
    ax.legend(loc='lower right', frameon=True)
    ax.set_ylim(0, 115)

    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            ax.annotate(f'{height}%',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontweight='bold', fontsize=9)

    autolabel(rects1)
    autolabel(rects2)

    plt.tight_layout()
    plt.savefig('Fig3_VCF_Benchmark_Performance.pdf', format='pdf', bbox_inches='tight')
    plt.savefig('Fig3_VCF_Benchmark_Performance.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()
    print("[✓] Generated Fig3_VCF_Benchmark_Performance")

if __name__ == '__main__':
    generate_figure_1()
    generate_figure_2()
    generate_figure_3()
    print("\nAll 3 figures successfully created in current workspace directory!")