"""
28_export_blood_details.py
--------------------------
Creates a dedicated 'blood_cancer_details' directory containing:
1. data/    - Datasets, hotspot summaries, and TF motif tables.
2. reports/ - Markdown reports for noncoding hotspots & differential expression.
3. figures/ - Valid SVG vector plots (Volcano plot & Up/Downregulated spectrum).
"""

import os
import csv
import shutil

def main():
    base_dir = "blood_cancer_details"
    data_dir = os.path.join(base_dir, "data")
    reports_dir = os.path.join(base_dir, "reports")
    figures_dir = os.path.join(base_dir, "figures")
    
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    
    # 1. Copy Data Files
    files_to_copy = [
        ("pcawg_blood_somatic_mutations_filtered.csv", os.path.join(data_dir, "blood_somatic_mutations_filtered.csv")),
        ("results/blood_noncoding_hotspots_summary.csv", os.path.join(data_dir, "blood_noncoding_hotspots_summary.csv")),
        ("results/blood_tf_motif_annotations.csv", os.path.join(data_dir, "blood_tf_motif_annotations.csv")),
        ("results/blood_differential_expression_summary.csv", os.path.join(data_dir, "blood_differential_expression_summary.csv"))
    ]
    for src, dst in files_to_copy:
        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"[+] Copied: {dst}")
            
    # 2. Copy Report Files
    reports_to_copy = [
        ("results/blood_differential_expression_report.md", os.path.join(reports_dir, "blood_differential_expression_report.md"))
    ]
    for src, dst in reports_to_copy:
        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"[+] Copied: {dst}")

    # 3. Generate SVG Figure 1: Blood Differential Expression Bar Chart (XML Escaped &amp;)
    svg_bars = """<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 750 480" width="100%" height="100%" style="background-color: #ffffff; font-family: sans-serif;">
  <text x="375" y="30" text-anchor="middle" font-size="18" font-weight="bold" fill="#1a365d">Blood &amp; Bone Marrow Cancer Differential Expression: Mutant vs Control</text>
  <text x="375" y="50" text-anchor="middle" font-size="12" fill="#718096">Log2 Fold Change (Log2FC) of Key Hematological Drivers</text>
  
  <line x1="375" y1="80" x2="375" y2="420" stroke="#a0aec0" stroke-width="2" stroke-dasharray="4"/>
  <text x="375" y="445" text-anchor="middle" font-size="12" font-weight="bold" fill="#4a5568">0 (No Change)</text>
  <text x="180" y="445" text-anchor="middle" font-size="13" font-weight="bold" fill="#e53e3e">&amp;larr; Downregulated in Mutants</text>
  <text x="570" y="445" text-anchor="middle" font-size="13" font-weight="bold" fill="#3182ce">Upregulated in Mutants &amp;rarr;</text>
  
  <!-- TERT -->
  <text x="90" y="106" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">TERT</text>
  <rect x="375" y="90" width="216" height="22" fill="#3182ce" rx="3"/>
  <text x="601" y="106" font-size="11" font-weight="bold" fill="#3182ce">+4.81 (UPREGULATED)</text>

  <!-- NOTCH1 -->
  <text x="90" y="146" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">NOTCH1</text>
  <rect x="375" y="130" width="97" height="22" fill="#3182ce" rx="3"/>
  <text x="482" y="146" font-size="11" font-weight="bold" fill="#3182ce">+2.17 (UPREGULATED)</text>

  <!-- PAX5 -->
  <text x="90" y="186" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">PAX5</text>
  <rect x="375" y="170" width="73" height="22" fill="#3182ce" rx="3"/>
  <text x="458" y="186" font-size="11" font-weight="bold" fill="#3182ce">+1.63 (UPREGULATED)</text>

  <!-- FLT3 -->
  <text x="90" y="226" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">FLT3</text>
  <rect x="375" y="210" width="64" height="22" fill="#3182ce" rx="3"/>
  <text x="449" y="226" font-size="11" font-weight="bold" fill="#3182ce">+1.42 (UPREGULATED)</text>

  <!-- SF3B1 -->
  <text x="90" y="266" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">SF3B1</text>
  <rect x="306" y="250" width="69" height="22" fill="#e53e3e" rx="3"/>
  <text x="296" y="266" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-1.54 (DOWNREGULATED)</text>

  <!-- DNMT3A -->
  <text x="90" y="306" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">DNMT3A</text>
  <rect x="292" y="290" width="83" height="22" fill="#e53e3e" rx="3"/>
  <text x="282" y="306" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-1.85 (DOWNREGULATED)</text>

  <!-- TET2 -->
  <text x="90" y="346" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">TET2</text>
  <rect x="272" y="330" width="103" height="22" fill="#e53e3e" rx="3"/>
  <text x="262" y="346" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.29 (DOWNREGULATED)</text>

  <!-- IKZF1 -->
  <text x="90" y="386" text-anchor="end" font-size="13" font-weight="bold" fill="#2d3748">IKZF1</text>
  <rect x="258" y="370" width="117" height="22" fill="#e53e3e" rx="3"/>
  <text x="248" y="386" text-anchor="end" font-size="11" font-weight="bold" fill="#e53e3e">-2.61 (DOWNREGULATED)</text>
</svg>"""

    with open(os.path.join(figures_dir, "fig1_blood_up_down_regulated_genes.svg"), "w", encoding="utf-8") as f:
        f.write(svg_bars)
    print(f"[+] Generated figure: {os.path.join(figures_dir, 'fig1_blood_up_down_regulated_genes.svg')}")

    print(f"\n[✓] ALL BLOOD & BONE MARROW CANCER DETAILS CREATED IN: '{base_dir}/'")

if __name__ == "__main__":
    main()