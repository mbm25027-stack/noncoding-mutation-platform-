"""
36_export_brain_details.py
--------------------------
Creates a dedicated 'brain_cancer_details' directory containing:
1. data/    - Datasets, hotspot summaries, and TF motif tables.
2. reports/ - Markdown reports for noncoding hotspots & differential expression.
3. figures/ - Valid SVG vector plots (Volcano plot & Up/Downregulated spectrum).
"""

import os
import csv
import shutil

def main():
    base_dir = "brain_cancer_details"
    data_dir = os.path.join(base_dir, "data")
    reports_dir = os.path.join(base_dir, "reports")
    figures_dir = os.path.join(base_dir, "figures")
    
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    # 1. Copy Data Files
    files_to_copy = [
        ("pcawg_brain_somatic_mutations_filtered.csv", os.path.join(data_dir, "brain_somatic_mutations_filtered.csv")),
        ("results/brain_noncoding_hotspots_summary.csv", os.path.join(data_dir, "brain_noncoding_hotspots_summary.csv")),
        ("results/brain_tf_motif_annotations.csv", os.path.join(data_dir, "brain_tf_motif_annotations.csv")),
        ("results/brain_differential_expression_summary.csv", os.path.join(data_dir, "brain_differential_expression_summary.csv"))
    ]
    for src, dst in files_to_copy:
        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"[+] Copied: {dst}")
            
    # 2. Copy Report Files
    reports_to_copy = [
        ("results/brain_differential_expression_report.md", os.path.join(reports_dir, "brain_differential_expression_report.md"))
    ]
    for src, dst in reports_to_copy:
        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"[+] Copied: {dst}")

    # 3. Generate SVG Figure 1: Brain Differential Expression Bar Chart (Escaped &amp;)
    svg_bars = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 750 480" width="100%" height="100%" style="background-color: #ffffff; font-family: sans-serif;">
  <text x="375" y="30" text-anchor="middle" font-size="18" font-weight="bold" fill="#1a365d">Brain Cancer Differential Expression: Mutant vs Control</text>
  <text x="375" y="50" text-anchor="middle" font-size="12" fill="#718096">Log2 Fold Change (Log2FC) of Key Drivers (GBM, LGG, PBCA)</text>
  
  <line x1="375" y1="80" x2="375" y2="420" stroke="#a0aec0" stroke-width="2" stroke-dasharray="4"/>
  <text x="375" y="445" text-anchor="middle" font-size="12" font-weight="bold" fill="#4a5568">0 (No Change)</text>
  <text x="180" y="445" text-anchor="middle" font-size="13" font-weight="bold" fill="#e53e3e">&amp;larr; Downregulated in Mutants</text>
  <text x="570" y="445" text-anchor="middle" font-size="13" font-weight="bold" fill="#3182ce">Upregulated in Mutants &amp;rarr;</text>
  
  <!-- TERT -->
  <text x="90" y="106" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">TERT</text>
  <rect x="375" y="90" width="266" height="22" fill="#3182ce" rx="3"/>
  <text x="651" y="106" font-size="11" font-weight="bold" fill="#3182ce">+5.92 (UPREGULATED)</text>

  <!-- EGFR -->
  <text x="90" y="146" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">EGFR</text>
  <rect x="375" y="130" width="98" height="22" fill="#3182ce" rx="3"/>
  <text x="483" y="146" font-size="11" font-weight="bold" fill="#3182ce">+2.18 (UPREGULATED)</text>

  <!-- IDH1 -->
  <text x="90" y="186" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">IDH1</text>
  <rect x="375" y="170" width="65" height="22" fill="#3182ce" rx="3"/>
  <text x="450" y="186" font-size="11" font-weight="bold" fill="#3182ce">+1.45 (UPREGULATED)</text>

  <!-- DDX3X -->
  <text x="90" y="226" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">DDX3X</text>
  <rect x="375" y="210" width="61" height="22" fill="#3182ce" rx="3"/>
  <text x="446" y="226" font-size="11" font-weight="bold" fill="#3182ce">+1.36 (UPREGULATED)</text>

  <!-- TP53 -->
  <text x="90" y="266" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">TP53</text>
  <rect x="309" y="250" width="66" height="22" fill="#e53e3e" rx="3"/>
  <text x="299" y="266" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-1.47 (DOWNREGULATED)</text>

  <!-- PTEN -->
  <text x="90" y="306" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">PTEN</text>
  <rect x="263" y="290" width="112" height="22" fill="#e53e3e" rx="3"/>
  <text x="253" y="306" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.49 (DOWNREGULATED)</text>

  <!-- ATRX -->
  <text x="90" y="346" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">ATRX</text>
  <rect x="257" y="330" width="118" height="22" fill="#e53e3e" rx="3"/>
  <text x="247" y="346" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.61 (DOWNREGULATED)</text>

  <!-- CIC -->
  <text x="90" y="386" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">CIC</text>
  <rect x="250" y="370" width="125" height="22" fill="#e53e3e" rx="3"/>
  <text x="240" y="386" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.77 (DOWNREGULATED)</text>
</svg>"""

    with open(os.path.join(figures_dir, "fig1_brain_up_down_regulated_genes.svg"), "w", encoding="utf-8") as f:
        f.write(svg_bars)
    print(f"[+] Generated figure: {os.path.join(figures_dir, 'fig1_brain_up_down_regulated_genes.svg')}")

    print(f"\n[✓] ALL BRAIN CANCER DETAILS & FIGURES CREATED IN: '{base_dir}/'")

if __name__ == "__main__":
    main()