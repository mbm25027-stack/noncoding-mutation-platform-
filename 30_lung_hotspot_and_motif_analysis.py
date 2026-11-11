"""
30_lung_hotspot_and_motif_analysis.py
-------------------------------------
Pure Standard Library script to identify recurrent noncoding hotspots
and transcription factor (TF) motif alterations in Lung Cancer (LUAD & LUSC).
"""

import csv
import os
from collections import defaultdict

def load_data(input_csv):
    with open(input_csv, mode='r', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def process_noncoding_hotspots(data):
    noncoding_categories = {"5'Flank", "5'UTR", "3'UTR", "Intron", "RNA", "Splice_Site"}
    nc_data = [row for row in data if row.get('Variant_Classification') in noncoding_categories]
    
    grouped = defaultdict(list)
    for row in nc_data:
        key = (
            row['Hugo_Symbol'],
            row['Chromosome'],
            row['Start_Position'],
            row['Variant_Classification'],
            row['HGVSp_Short']
        )
        grouped[key].append(row)
    
    hotspots = []
    for key, rows in grouped.items():
        gene, chrom, pos, vclass, hgvsp = key
        donors = set(r['Donor_ID'] for r in rows if r.get('Donor_ID'))
        projects = set(r['Project_Code'] for r in rows if r.get('Project_Code'))
        vafs = [float(r['VAF']) for r in rows if r.get('VAF') not in (None, '', 'N/A')]
        mean_vaf = round(sum(vafs) / len(vafs), 3) if vafs else 0.0
        
        hotspots.append({
            'Hugo_Symbol': gene,
            'Chromosome': chrom,
            'Start_Position': pos,
            'Variant_Classification': vclass,
            'HGVSp_Short': hgvsp,
            'Donor_Count': len(donors),
            'Projects': ', '.join(projects),
            'Mean_VAF': mean_vaf
        })
    
    hotspots.sort(key=lambda x: x['Donor_Count'], reverse=True)
    return nc_data, hotspots

def annotate_tf_motifs(nc_data):
    motif_annotations = []
    for row in nc_data:
        gene = row.get('Hugo_Symbol', '')
        hgvsp = str(row.get('HGVSp_Short', ''))
        vclass = row.get('Variant_Classification', '')
        pos = str(row.get('Start_Position', ''))
        
        tf_effect = "Unknown / Neutral"
        predicted_tf = "None"
        motif_change = "N/A"
        
        if gene == "TERT" or "TERT_p" in hgvsp:
            if "C228T" in hgvsp or "1295228" in pos:
                tf_effect = "De novo motif creation (GABP / ETS binding)"
                predicted_tf = "GABPA/ETS1"
                motif_change = "CCCCTT -> CCGGAA"
            elif "C250T" in hgvsp or "1295250" in pos:
                tf_effect = "De novo motif creation (GABP / ETS binding)"
                predicted_tf = "GABPA/ETS1"
                motif_change = "CCCCTT -> CCGGAA"
        elif gene == "EGFR" and vclass == "5'Flank":
            tf_effect = "Promoter Sp1/E2F binding enhancement"
            predicted_tf = "Sp1 / E2F1"
            motif_change = "GC-box alteration"
        elif gene == "NKX2-1" and vclass == "5'Flank":
            tf_effect = "Lineage transcription factor binding gain"
            predicted_tf = "FOXA1 / TTF-1"
            motif_change = "Enhancer motif shift"
        elif gene == "MET" and vclass == "Splice_Site":
            tf_effect = "Splice site donor disruption (Exon 14 skipping)"
            predicted_tf = "U1/U2 snRNP"
            motif_change = "Exon 14 splice loss"
            
        vaf_val = float(row['VAF']) if row.get('VAF') not in (None, '', 'N/A') else 0.0
        
        motif_annotations.append({
            'Donor_ID': row.get('Donor_ID', ''),
            'Project_Code': row.get('Project_Code', ''),
            'Gene': gene,
            'Chromosome': row.get('Chromosome', ''),
            'Position': pos,
            'Variant_Classification': vclass,
            'HGVSp_Short': hgvsp,
            'VAF': vaf_val,
            'TF_Effect': tf_effect,
            'Predicted_TF': predicted_tf,
            'Motif_Change': motif_change
        })
    return motif_annotations

def save_csv(data, filename, fieldnames):
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def main():
    input_file = "pcawg_lung_somatic_mutations_filtered.csv"
    output_hotspots = "results/lung_noncoding_hotspots_summary.csv"
    output_tf = "results/lung_tf_motif_annotations.csv"
    
    if not os.path.exists(input_file):
        print(f"[!] File not found: {input_file}")
        return
        
    data = load_data(input_file)
    nc_data, hotspots = process_noncoding_hotspots(data)
    tf_data = annotate_tf_motifs(nc_data)
    
    hotspot_fields = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Variant_Classification', 'HGVSp_Short', 'Donor_Count', 'Projects', 'Mean_VAF']
    tf_fields = ['Donor_ID', 'Project_Code', 'Gene', 'Chromosome', 'Position', 'Variant_Classification', 'HGVSp_Short', 'VAF', 'TF_Effect', 'Predicted_TF', 'Motif_Change']
    
    save_csv(hotspots, output_hotspots, hotspot_fields)
    save_csv(tf_data, output_tf, tf_fields)
    
    print(f"[+] Saved lung noncoding hotspot summary to: {output_hotspots}")
    print(f"[+] Saved lung TF motif analysis to: {output_tf}")

if __name__ == "__main__":
    main()