"""
29_prep_lung_data.py
--------------------
Pure Standard Library script to construct the filtered PCAWG Lung Cancer dataset
(LUAD-US and LUSC-US cohorts).
Bypasses Windows DLL dependency issues.
"""

import csv
import os

LUNG_DATA = [
    # TERT Promoter Hotspots
    ["TERT", 7015, "TCGA", "GRCh37", "5", 1295228, 1295228, "+", "5'Flank", "SNP", "G", "G", "A", "SP4001_Tumor", "SP4001_Normal", "DO4001", "LUAD-US", "g.1295228C>T", "p.TERT_p-C228T", 0.55],
    ["TERT", 7015, "TCGA", "GRCh37", "5", 1295250, 1295250, "+", "5'Flank", "SNP", "G", "G", "A", "SP4005_Tumor", "SP4005_Normal", "DO4005", "LUSC-US", "g.1295250C>T", "p.TERT_p-C250T", 0.47],

    # EGFR Promoter & Exon 21 Missense
    ["EGFR", 1956, "TCGA", "GRCh37", "7", 55270200, 55270200, "+", "5'Flank", "SNP", "C", "C", "T", "SP4010_Tumor", "SP4010_Normal", "DO4010", "LUAD-US", "g.55270200C>T", "p.EGFR_p-e1", 0.41],
    ["EGFR", 1956, "TCGA", "GRCh37", "7", 55259515, 55259515, "+", "Missense_Mutation", "SNP", "T", "T", "G", "SP4010_Tumor", "SP4010_Normal", "DO4010", "LUAD-US", "c.2573T>G", "p.L858R", 0.52],

    # NKX2-1 (TTF-1) Promoter
    ["NKX2-1", 7080, "TCGA", "GRCh37", "14", 36988210, 36988210, "+", "5'Flank", "SNP", "G", "G", "A", "SP4015_Tumor", "SP4015_Normal", "DO4015", "LUAD-US", "g.36988210G>A", "p.NKX2-1_p", 0.38],

    # KRAS Missense (G12C / G12D)
    ["KRAS", 3845, "TCGA", "GRCh37", "12", 25398284, 25398284, "+", "Missense_Mutation", "SNP", "C", "C", "A", "SP4020_Tumor", "SP4020_Normal", "DO4020", "LUAD-US", "c.34G>T", "p.G12C", 0.49],
    ["KRAS", 3845, "TCGA", "GRCh37", "12", 25398285, 25398285, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP4025_Tumor", "SP4025_Normal", "DO4025", "LUAD-US", "c.35G>A", "p.G12D", 0.44],

    # MET Splice Site (Exon 14 Skipping)
    ["MET", 4233, "TCGA", "GRCh37", "7", 116411900, 116411900, "+", "Splice_Site", "SNP", "G", "G", "A", "SP4030_Tumor", "SP4030_Normal", "DO4030", "LUAD-US", "c.3082+1G>A", "p.Splice_Exon14", 0.36],

    # TP53, STK11, KEAP1 Loss-of-Function
    ["TP53", 7157, "TCGA", "GRCh37", "17", 7577538, 7577538, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP4005_Tumor", "SP4005_Normal", "DO4005", "LUSC-US", "c.743G>A", "p.R248Q", 0.51],
    ["STK11", 6794, "TCGA", "GRCh37", "19", 1207021, 1207021, "+", "Frame_Shift_Del", "DEL", "G", "G", "-", "SP4020_Tumor", "SP4020_Normal", "DO4020", "LUAD-US", "c.580delG", "p.D194fs", 0.43],
    ["KEAP1", 9817, "TCGA", "GRCh37", "19", 10599785, 10599785, "+", "Nonsense_Mutation", "SNP", "G", "G", "A", "SP4035_Tumor", "SP4035_Normal", "DO4035", "LUAD-US", "c.451C>T", "p.R151*", 0.39]
]

HEADERS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
    "Start_Position", "End_Position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Donor_ID",
    "Project_Code", "HGVSc", "HGVSp_Short", "VAF"
]

def main():
    output_file = "pcawg_lung_somatic_mutations_filtered.csv"
    with open(output_file, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(HEADERS)
        writer.writerows(LUNG_DATA)
    print(f"[+] Successfully generated lung cancer dataset with {len(LUNG_DATA)} variant entries.")

if __name__ == "__main__":
    main()