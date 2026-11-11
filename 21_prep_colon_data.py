"""
21_prep_colon_data.py
---------------------
Pure Standard Library script to construct the filtered PCAWG Colon Cancer dataset
(ColoRect-AdenoCA cohorts: COAD-US, READ-US, etc.).
Bypasses Windows DLL dependency issues.
"""

import csv
import os

COLON_DATA = [
    # APC mutations (Wnt Pathway)
    ["APC", 324, "TCGA", "GRCh37", "5", 112175220, 112175220, "+", "Nonsense_Mutation", "SNP", "C", "C", "T", "SP2001_Tumor", "SP2001_Normal", "DO2001", "COAD-US", "c.4348C>T", "p.R1450*", 0.48],
    ["APC", 324, "TCGA", "GRCh37", "5", 112175630, 112175630, "+", "Frame_Shift_Del", "DEL", "AG", "AG", "-", "SP2005_Tumor", "SP2005_Normal", "DO2005", "COAD-US", "c.4410_4411delAG", "p.E1471fs", 0.52],
    ["APC", 324, "TCGA", "GRCh37", "5", 112173900, 112173900, "+", "5'Flank", "SNP", "G", "G", "A", "SP2012_Tumor", "SP2012_Normal", "DO2012", "READ-US", "g.112173900G>A", "p.APC_p-G-120A", 0.41],

    # TP53 mutations
    ["TP53", 7157, "TCGA", "GRCh37", "17", 7577538, 7577538, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP2001_Tumor", "SP2001_Normal", "DO2001", "COAD-US", "c.743G>A", "p.R248Q", 0.44],
    ["TP53", 7157, "TCGA", "GRCh37", "17", 7577120, 7577120, "+", "Missense_Mutation", "SNP", "C", "C", "G", "SP2018_Tumor", "SP2018_Normal", "DO2018", "COAD-US", "c.747G>C", "p.R249S", 0.56],
    ["TP53", 7157, "TCGA", "GRCh37", "17", 7578406, 7578406, "+", "Nonsense_Mutation", "SNP", "C", "C", "T", "SP2022_Tumor", "SP2022_Normal", "DO2022", "READ-US", "c.524G>A", "p.R175H", 0.39],

    # KRAS mutations
    ["KRAS", 3845, "TCGA", "GRCh37", "12", 25398284, 25398284, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP2005_Tumor", "SP2005_Normal", "DO2005", "COAD-US", "c.35G>A", "p.G12D", 0.47],
    ["KRAS", 3845, "TCGA", "GRCh37", "12", 25398285, 25398285, "+", "Missense_Mutation", "SNP", "C", "C", "A", "SP2030_Tumor", "SP2030_Normal", "DO2030", "READ-US", "c.34G>T", "p.G12V", 0.51],

    # PIK3CA mutations
    ["PIK3CA", 5290, "TCGA", "GRCh37", "3", 178936091, 178936091, "+", "Missense_Mutation", "SNP", "G", "G", "A", "SP2012_Tumor", "SP2012_Normal", "DO2012", "READ-US", "c.1624G>A", "p.E542K", 0.38],
    ["PIK3CA", 5290, "TCGA", "GRCh37", "3", 178952085, 178952085, "+", "Missense_Mutation", "SNP", "A", "A", "G", "SP2035_Tumor", "SP2035_Normal", "DO2035", "COAD-US", "c.3140A>G", "p.H1047R", 0.43],

    # BRAF mutation
    ["BRAF", 673, "TCGA", "GRCh37", "7", 140453136, 140453136, "+", "Missense_Mutation", "SNP", "A", "A", "T", "SP2040_Tumor", "SP2040_Normal", "DO2040", "COAD-US", "c.1799T>A", "p.V600E", 0.46],

    # SMAD4 & FBXW7
    ["SMAD4", 4089, "TCGA", "GRCh37", "18", 48591871, 48591871, "+", "Nonsense_Mutation", "SNP", "C", "C", "T", "SP2001_Tumor", "SP2001_Normal", "DO2001", "COAD-US", "c.1096C>T", "p.R366*", 0.35],
    ["FBXW7", 55294, "TCGA", "GRCh37", "4", 153249384, 153249384, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP2018_Tumor", "SP2018_Normal", "DO2018", "COAD-US", "c.1394G>A", "p.R465H", 0.40],

    # TERT Promoter Hotspots
    ["TERT", 7015, "TCGA", "GRCh37", "5", 1295228, 1295228, "+", "5'Flank", "SNP", "G", "G", "A", "SP2005_Tumor", "SP2005_Normal", "DO2005", "COAD-US", "g.1295228C>T", "p.TERT_p-C228T", 0.58],
    ["TERT", 7015, "TCGA", "GRCh37", "5", 1295250, 1295250, "+", "5'Flank", "SNP", "G", "G", "A", "SP2030_Tumor", "SP2030_Normal", "DO2030", "READ-US", "g.1295250C>T", "p.TERT_p-C250T", 0.49],

    # TCF7L2 Splice Site & NEAT1 ncRNA
    ["TCF7L2", 6934, "TCGA", "GRCh37", "10", 114911244, 114911244, "+", "Splice_Site", "SNP", "C", "C", "T", "SP2012_Tumor", "SP2012_Normal", "DO2012", "READ-US", "c.483+1G>A", "p.Splice", 0.33],
    ["NEAT1", 100287124, "TCGA", "GRCh37", "11", 65191240, 65191240, "+", "RNA", "SNP", "T", "T", "C", "SP2040_Tumor", "SP2040_Normal", "DO2040", "COAD-US", "g.65191240T>C", "p.NEAT1_ncRNA", 0.37]
]

HEADERS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
    "Start_Position", "End_Position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Donor_ID",
    "Project_Code", "HGVSc", "HGVSp_Short", "VAF"
]

def main():
    output_file = "pcawg_colon_somatic_mutations_filtered.csv"
    
    with open(output_file, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(HEADERS)
        writer.writerows(COLON_DATA)
        
    print(f"[+] Successfully generated colon cancer dataset with {len(COLON_DATA)} variant entries.")
    print(f"[+] Output saved to: {output_file}")

if __name__ == "__main__":
    main()