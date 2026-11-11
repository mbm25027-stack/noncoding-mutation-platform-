"""
17_plot_liver_summary.py
-----------------------
Pure Standard Library HTML/SVG visualizer for PCAWG Liver Cancer Noncoding Hotspots.
Generates an interactive visual dashboard in results/liver_visual_summary.html.
"""

import csv
import os

def load_csv(filepath):
    if not os.path.exists(filepath):
        return []
    with open(filepath, mode='r', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def main():
    hotspots = load_csv("results/liver_noncoding_hotspots_summary.csv")
    tf_data = load_csv("results/liver_tf_motif_annotations.csv")
    
    html_content = """<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>PCAWG Liver Noncoding Analysis Dashboard</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; margin: 30px; background: #f8f9fa; color: #333; }
        h1 { color: #1a365d; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; }
        h2 { color: #2b6cb0; margin-top: 30px; }
        .card { background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); margin-bottom: 25px; }
        table { width: 100%; border-collapse: collapse; margin-top: 10px; }
        th, td { text-align: left; padding: 12px; border-bottom: 1px solid #e2e8f0; }
        th { background-color: #ebf8ff; color: #2c5282; font-weight: 600; }
        tr:hover { background-color: #f7fafc; }
        .badge-hotspot { background: #e53e3e; color: white; padding: 3px 8px; border-radius: 4px; font-size: 0.85em; }
        .badge-tf { background: #3182ce; color: white; padding: 3px 8px; border-radius: 4px; font-size: 0.85em; }
        .bar-container { background: #edf2f7; width: 100%; height: 20px; border-radius: 10px; overflow: hidden; }
        .bar-fill { background: #4299e1; height: 100%; text-align: right; padding-right: 5px; color: white; font-size: 12px; line-height: 20px; }
    </style>
</head>
<body>
    <h1>PCAWG Liver Cancer Noncoding Analysis Dashboard</h1>
    
    <div class="card">
        <h2>1. Recurrent Noncoding Promoter & Splice Hotspots</h2>
        <table>
            <thead>
                <tr>
                    <th>Hugo Symbol</th>
                    <th>Genomic Position</th>
                    <th>Classification</th>
                    <th>Recurrent Donors</th>
                    <th>Mean VAF</th>
                </tr>
            </thead>
            <tbody>
"""
    for h in hotspots:
        vaf_pct = float(h['Mean_VAF']) * 100
        html_content += f"""
                <tr>
                    <td><strong>{h['Hugo_Symbol']}</strong></td>
                    <td>Chr{h['Chromosome']}:{h['Start_Position']}</td>
                    <td><span class="badge-hotspot">{h['Variant_Classification']}</span></td>
                    <td>{h['Donor_Count']}</td>
                    <td>
                        <div class="bar-container">
                            <div class="bar-fill" style="width: {min(vaf_pct*1.2, 100)}%;">{h['Mean_VAF']}</div>
                        </div>
                    </td>
                </tr>"""

    html_content += """
            </tbody>
        </table>
    </div>

    <div class="card">
        <h2>2. Transcription Factor (TF) Motif Alteration Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Gene</th>
                    <th>HGVS Short</th>
                    <th>Predicted TF</th>
                    <th>Functional Effect</th>
                    <th>Motif Sequence Shift</th>
                </tr>
            </thead>
            <tbody>
"""
    for t in tf_data:
        html_content += f"""
                <tr>
                    <td><strong>{t['Gene']}</strong></td>
                    <td>{t['HGVSp_Short']}</td>
                    <td><span class="badge-tf">{t['Predicted_TF']}</span></td>
                    <td>{t['TF_Effect']}</td>
                    <td><code>{t['Motif_Change']}</code></td>
                </tr>"""

    html_content += """
            </tbody>
        </table>
    </div>
</body>
</html>
"""

    output_path = "results/liver_visual_summary.html"
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    
    print(f"[+] Interactive visual summary created at: {output_path}")

if __name__ == "__main__":
    main()