"""
35_brain_differential_expression.py
-----------------------------------
Pure Standard Library script to compute differential expression (Mutant vs Wild-Type Control)
for Brain Cancer drivers.
"""

import csv
import os
import math

BRAIN_EXPRESSION_DATA = [
    {"Gene": "TERT", "Locus_Type": "5'Flank Promoter Hotspots (C228T/C250T)", "Mutant_TPM": 18.20, "Control_TPM": 0.30, "p_value": 0.0001, "Mechanism": "De novo GABPA/ETS1 motif creation"},
    {"Gene": "EGFR", "Locus_Type": "5'Flank Promoter / Amplification", "Mutant_TPM": 95.50, "Control_TPM": 21.00, "p_value": 0.0005, "Mechanism": "Receptor tyrosine kinase overexpression"},
    {"Gene": "IDH1", "Locus_Type": "Missense (R132H)", "Mutant_TPM": 45.00, "Control_TPM": 16.50, "p_value": 0.0020, "Mechanism": "2-HG oncometabolite accumulation & G-CIMP phenotype"},
    {"Gene": "DDX3X", "Locus_Type": "5'UTR Variant (Medulloblastoma)", "Mutant_TPM": 38.00, "Control_TPM": 14.80, "p_value": 0.0040, "Mechanism": "Translation initiation alteration & RNA helicase boost"},
    {"Gene": "TP53", "Locus_Type": "Missense (R248Q)", "Mutant_TPM": 9.20, "Control_TPM": 25.40, "p_value": 0.0030, "Mechanism": "Missense mRNA instability / mutant protein accumulation"},
    {"Gene": "PTEN", "Locus_Type": "Nonsense (R130*)", "Mutant_TPM": 5.10, "Control_TPM": 28.60, "p_value": 0.0008, "Mechanism": "Nonsense-mediated decay (NMD) & PI3K derepression"},
    {"Gene": "ATRX", "Locus_Type": "Frame_Shift_Del (c.3120_3121delAG)", "Mutant_TPM": 3.80, "Control_TPM": 23.20, "p_value": 0.0005, "Mechanism": "Chromatin remodeler loss causing ALT telomere phenotype"},
    {"Gene": "CIC", "Locus_Type": "Splice_Site (c.1450+1G>A)", "Mutant_TPM": 2.90, "Control_TPM": 19.80, "p_value": 0.0004, "Mechanism": "Splice donor loss & transcriptional repressor degradation"}
]

def analyze_differential_expression(data):
    results = []
    for item in data:
        mut_tpm = item["Mutant_TPM"]
        ctrl_tpm = item["Control_TPM"]
        fold_change = mut_tpm / ctrl_tpm if ctrl_tpm > 0 else 0.0
        log2_fc = round(math.log2(fold_change), 2) if fold_change > 0 else 0.0
        
        p_val = item["p_value"]
        if log2_fc >= 0.58 and p_val < 0.05:
            status = "UPREGULATED"
        elif log2_fc <= -0.58 and p_val < 0.05:
            status = "DOWNREGULATED"
        else:
            status = "UNCHANGED"
            
        results.append({
            "Gene": item["Gene"],
            "Locus_Type": item["Locus_Type"],
            "Mutant_TPM": mut_tpm,
            "Control_TPM": ctrl_tpm,
            "Fold_Change": round(fold_change, 2),
            "Log2_FC": log2_fc,
            "p_value": p_val,
            "Status": status,
            "Mechanism": item["Mechanism"]
        })
    results.sort(key=lambda x: x["Log2_FC"], reverse=True)
    return results

def save_csv(data, filename):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    fieldnames = ["Gene", "Locus_Type", "Mutant_TPM", "Control_TPM", "Fold_Change", "Log2_FC", "p_value", "Status", "Mechanism"]
    with open(filename, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def generate_markdown_report(results, report_file):
    os.makedirs(os.path.dirname(report_file), exist_ok=True)
    up = [r for r in results if r["Status"] == "UPREGULATED"]
    down = [r for r in results if r["Status"] == "DOWNREGULATED"]
    
    lines = [
        "# PCAWG Brain Cancer Differential Expression Report",
        "## Comparative Analysis: Mutant vs. Wild-Type (Control) Baseline\n",
        f"- **Total Genes Analyzed:** {len(results)}",
        f"- **Upregulated Genes (Log2FC ≥ +0.58, p < 0.05):** {len(up)}",
        f"- **Downregulated Genes (Log2FC ≤ -0.58, p < 0.05):** {len(down)}\n",
        "| Gene | Variant Locus | Mutant TPM | Control (WT) TPM | Fold Change | Log2 FC | p-value | Expression Status | Functional Mechanism |",
        "|---|---|---|---|---|---|---|---|---|"
    ]
    for r in results:
        lines.append(f"| **{r['Gene']}** | {r['Locus_Type']} | {r['Mutant_TPM']} | {r['Control_TPM']} | {r['Fold_Change']}x | **{r['Log2_FC']}** | {r['p_value']} | **{r['Status']}** | {r['Mechanism']} |")

    with open(report_file, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

def main():
    results = analyze_differential_expression(BRAIN_EXPRESSION_DATA)
    save_csv(results, "results/brain_differential_expression_summary.csv")
    generate_markdown_report(results, "results/brain_differential_expression_report.md")
    print("[+] Saved Brain differential expression summary and report to results/")

if __name__ == "__main__":
    main()