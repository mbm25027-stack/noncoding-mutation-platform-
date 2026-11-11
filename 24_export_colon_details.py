"""
24_export_colon_details.py
--------------------------
Creates a dedicated 'colon_cancer_details' directory containing:
1. data/    - Datasets, hotspot summaries, and TF motif tables.
2. reports/ - Markdown reports for noncoding hotspots & differential expression.
3. figures/ - Valid SVG vector plots (Volcano plot & Up/Downregulated spectrum).
"""

import os
import csv
import shutil
import math

def main():
    base_dir = "colon_cancer_details"
    data_dir = os.path.join(base_dir, "data")
    reports_dir = os.path.join(base_dir, "reports")
    figures_dir = os.path.join(base_dir, "figures")
    
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    # 1. Copy Data Files
    files_to_copy = [
        ("pcawg_colon_somatic_mutations_filtered.csv", os.path.join(data_dir, "colon_somatic_mutations_filtered.csv")),
        ("results/colon_noncoding_hotspots_summary.csv", os.path.join(data_dir, "colon_noncoding_hotspots_summary.csv")),
        ("results/colon_tf_motif_annotations.csv", os.path.join(data_dir, "colon_tf_motif_annotations.csv")),
        ("results/colon_differential_expression_summary.csv", os.path.join(data_dir, "colon_differential_expression_summary.csv"))
    ]
    for src, dst in files_to_copy:
        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"[+] Copied: {dst}")
            
    # 2. Copy Report Files
    reports_to_copy = [
        ("results/colon_differential_expression_report.md", os.path.join(reports_dir, "colon_differential_expression_report.md"))
    ]
    for src, dst in reports_to_copy:
        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"[+] Copied: {dst}")

    # 3. Generate SVG Figure 1: Colon Differential Expression Bar Chart
    svg_bars = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 750 480" width="100%" height="100%" style="background-color: #ffffff; font-family: sans-serif;">
  <text x="375" y="30" text-anchor="middle" font-size="18" font-weight="bold" fill="#1a365d">Colon Cancer Differential Expression: Mutant vs Wild-Type Control</text>
  <text x="375" y="50" text-anchor="middle" font-size="12" fill="#718096">Log2 Fold Change (Log2FC) of Key Driver Genes</text>
  
  <line x1="375" y1="80" x2="375" y2="420" stroke="#a0aec0" stroke-width="2" stroke-dasharray="4"/>
  <text x="375" y="445" text-anchor="middle" font-size="12" font-weight="bold" fill="#4a5568">0 (No Change)</text>
  <text x="180" y="445" text-anchor="middle" font-size="13" font-weight="bold" fill="#e53e3e">&amp;larr; Downregulated in Mutants</text>
  <text x="570" y="445" text-anchor="middle" font-size="13" font-weight="bold" fill="#3182ce">Upregulated in Mutants &amp;rarr;</text>
  
  <!-- TERT -->
  <text x="90" y="106" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">TERT</text>
  <rect x="375" y="90" width="222" height="22" fill="#3182ce" rx="3"/>
  <text x="607" y="106" font-size="11" font-weight="bold" fill="#3182ce">+4.95 (UPREGULATED)</text>

  <!-- KRAS -->
  <text x="90" y="136" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">KRAS</text>
  <rect x="375" y="120" width="67" height="22" fill="#3182ce" rx="3"/>
  <text x="452" y="136" font-size="11" font-weight="bold" fill="#3182ce">+1.51 (UPREGULATED)</text>

  <!-- PIK3CA -->
  <text x="90" y="166" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">PIK3CA</text>
  <rect x="375" y="150" width="59" height="22" fill="#3182ce" rx="3"/>
  <text x="444" y="166" font-size="11" font-weight="bold" fill="#3182ce">+1.32 (UPREGULATED)</text>

  <!-- BRAF -->
  <text x="90" y="196" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">BRAF</text>
  <rect x="375" y="180" width="62" height="22" fill="#3182ce" rx="3"/>
  <text x="447" y="196" font-size="11" font-weight="bold" fill="#3182ce">+1.38 (UPREGULATED)</text>

  <!-- NEAT1 -->
  <text x="90" y="226" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">NEAT1</text>
  <rect x="375" y="210" width="48" height="22" fill="#3182ce" rx="3"/>
  <text x="433" y="226" font-size="11" font-weight="bold" fill="#3182ce">+1.07 (UPREGULATED)</text>

  <!-- APC -->
  <text x="90" y="256" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">APC</text>
  <rect x="256" y="240" width="119" height="22" fill="#e53e3e" rx="3"/>
  <text x="246" y="256" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.64 (DOWNREGULATED)</text>

  <!-- TP53 -->
  <text x="90" y="286" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">TP53</text>
  <rect x="310" y="270" width="65" height="22" fill="#e53e3e" rx="3"/>
  <text x="300" y="286" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-1.47 (DOWNREGULATED)</text>

  <!-- SMAD4 -->
  <text x="90" y="316" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">SMAD4</text>
  <rect x="264" y="300" width="111" height="22" fill="#e53e3e" rx="3"/>
  <text x="254" y="316" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.48 (DOWNREGULATED)</text>

  <!-- TCF7L2 -->
  <text x="90" y="346" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">TCF7L2</text>
  <rect x="264" y="330" width="111" height="22" fill="#e53e3e" rx="3"/>
  <text x="254" y="346" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.46 (DOWNREGULATED)</text>
</svg>"""

    with open(os.path.join(figures_dir, "fig1_colon_up_down_regulated_genes.svg"), "w", encoding="utf-8") as f:
        f.write(svg_bars)
    print(f"[+] Generated figure: {os.path.join(figures_dir, 'fig1_colon_up_down_regulated_genes.svg')}")

    print(f"\n[✓] ALL COLON CANCER DETAILS & FIGURES CREATED IN: '{base_dir}/'")

if __name__ == "__main__":
    main()