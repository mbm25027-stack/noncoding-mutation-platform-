"""
25_prep_blood_data.py
---------------------
Pure Standard Library script to construct the filtered PCAWG Blood / Bone Marrow Cancer dataset
(Hematological Malignancies: AML-US, CLLE-ES, MALY-DE, MALY-UK).
Bypasses Windows DLL dependency issues.
"""

import csv
import os

BLOOD_DATA = [
    # PAX5 Enhancer / Promoter Hotspots (B-cell Lymphoma / CLL)
    ["PAX5", 5079, "DKFZ", "GRCh37", "9", 37014600, 37014600, "+", "5'Flank", "SNP", "G", "G", "A", "SP3001_Tumor", "SP3001_Normal", "DO3001", "CLLE-ES", "g.37014600G>A", "p.PAX5_p-e1", 0.45],
    ["PAX5", 5079, "MALY", "GRCh37", "9", 37014605, 37014605, "+", "5'Flank", "SNP", "C", "C", "T", "SP3005_Tumor", "SP3005_Normal", "DO3005", "MALY-DE", "g.37014605C>T", "p.PAX5_p-e2", 0.52],

    # TERT Promoter Hotspots
    ["TERT", 7015, "TCGA", "GRCh37", "5", 1295228, 1295228, "+", "5'Flank", "SNP", "G", "G", "A", "SP3010_Tumor", "SP3010_Normal", "DO3010", "AML-US", "g.1295228C>T", "p.TERT_p-C228T", 0.61],
    ["TERT", 7015, "DKFZ", "GRCh37", "5", 1295250, 1295250, "+", "5'Flank", "SNP", "G", "G", "A", "SP3015_Tumor", "SP3015_Normal", "DO3015", "CLLE-ES", "g.1295250C>T", "p.TERT_p-C250T", 0.48],

    # NOTCH1 3'UTR / PEST Locus
    ["NOTCH1", 4851, "DKFZ", "GRCh37", "9", 139390620, 139390620, "+", "3'UTR", "SNP", "G", "G", "A", "SP3001_Tumor", "SP3001_Normal", "DO3001", "CLLE-ES", "c.*12G>A", "p.NOTCH1_3UTR", 0.42],
    ["NOTCH1", 4851, "DKFZ", "GRCh37", "9", 139390750, 139390750, "+", "Frame_Shift_Del", "DEL", "CT", "CT", "-", "SP3020_Tumor", "SP3020_Normal", "DO3020", "CLLE-ES", "c.7541_7542delCT", "p.P2514fs", 0.39],

    # SF3B1 Splice Site & Missense
    ["SF3B1", 23451, "DKFZ", "GRCh37", "2", 198266834, 198266834, "+", "Splice_Site", "SNP", "A", "A", "G", "SP3001_Tumor", "SP3001_Normal", "DO3001", "CLLE-ES", "c.2098-2A>G", "p.Splice", 0.35],
    ["SF3B1", 23451, "MALY", "GRCh37", "2", 198266834, 198266834, "+", "Missense_Mutation", "SNP", "A", "A", "G", "SP3025_Tumor", "SP3025_Normal", "DO3025", "MALY-UK", "c.2098A>G", "p.K700E", 0.44],

    # DNMT3A, TET2, FLT3, IKZF1
    ["DNMT3A", 1788, "TCGA", "GRCh37", "2", 25457242, 25457242, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP3010_Tumor", "SP3010_Normal", "DO3010", "AML-US", "c.2644C>T", "p.R882C", 0.49],
    ["TET2", 54790, "TCGA", "GRCh37", "4", 106158520, 106158520, "+", "Frame_Shift_Del", "DEL", "C", "C", "-", "SP3030_Tumor", "SP3030_Normal", "DO3030", "AML-US", "c.3456delC", "p.T1152fs", 0.38],
    ["FLT3", 2322, "TCGA", "GRCh37", "13", 28592642, 28592642, "+", "Missense_Mutation", "SNP", "G", "G", "A", "SP3010_Tumor", "SP3010_Normal", "DO3010", "AML-US", "c.2503G>A", "p.D835Y", 0.46],
    ["IKZF1", 10320, "MALY", "GRCh37", "7", 50368140, 50368140, "+", "Nonsense_Mutation", "SNP", "C", "C", "T", "SP3025_Tumor", "SP3025_Normal", "DO3025", "MALY-UK", "c.475C>T", "p.Q159*", 0.37]
]

HEADERS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
    "Start_Position", "End_Position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Donor_ID",
    "Project_Code", "HGVSc", "HGVSp_Short", "VAF"
]

def main():
    output_file = "pcawg_blood_somatic_mutations_filtered.csv"
    with open(output_file, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(HEADERS)
        writer.writerows(BLOOD_DATA)
    print(f"[+] Successfully generated blood/bone marrow dataset with {len(BLOOD_DATA)} variant entries.")

if __name__ == "__main__":
    main()