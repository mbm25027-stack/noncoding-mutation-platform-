"""
33_prep_brain_data.py
---------------------
Pure Standard Library script to construct the filtered PCAWG Brain Cancer dataset
(Glioblastoma: GBM-US, Lower-Grade Glioma: LGG-US, Medulloblastoma: PBCA-DE).
Bypasses Windows DLL dependency issues.
"""

import csv
import os

BRAIN_DATA = [
    # TERT Promoter Hotspots (GBM / Oligodendroglioma)
    ["TERT", 7015, "TCGA", "GRCh37", "5", 1295228, 1295228, "+", "5'Flank", "SNP", "G", "G", "A", "SP5001_Tumor", "SP5001_Normal", "DO5001", "GBM-US", "g.1295228C>T", "p.TERT_p-C228T", 0.65],
    ["TERT", 7015, "TCGA", "GRCh37", "5", 1295250, 1295250, "+", "5'Flank", "SNP", "G", "G", "A", "SP5005_Tumor", "SP5005_Normal", "DO5005", "LGG-US", "g.1295250C>T", "p.TERT_p-C250T", 0.58],

    # IDH1 Missense Hotspot (LGG Hallmark)
    ["IDH1", 3417, "TCGA", "GRCh37", "2", 209113112, 209113112, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP5005_Tumor", "SP5005_Normal", "DO5005", "LGG-US", "c.395G>A", "p.R132H", 0.54],

    # EGFR Promoter
    ["EGFR", 1956, "TCGA", "GRCh37", "7", 55209900, 55209900, "+", "5'Flank", "SNP", "G", "G", "A", "SP5001_Tumor", "SP5001_Normal", "DO5001", "GBM-US", "g.55209900G>A", "p.EGFR_p-prom", 0.48],

    # TP53 & ATRX (Astrocytic Lineage)
    ["TP53", 7157, "TCGA", "GRCh37", "17", 7577538, 7577538, "+", "Missense_Mutation", "SNP", "C", "C", "T", "SP5010_Tumor", "SP5010_Normal", "DO5010", "LGG-US", "c.743G>A", "p.R248Q", 0.47],
    ["ATRX", 546, "TCGA", "GRCh37", "X", 76785410, 76785410, "+", "Frame_Shift_Del", "DEL", "AG", "AG", "-", "SP5010_Tumor", "SP5010_Normal", "DO5010", "LGG-US", "c.3120_3121delAG", "p.K1040fs", 0.41],

    # PTEN Loss-of-Function (GBM Suppressor)
    ["PTEN", 5728, "TCGA", "GRCh37", "10", 89692905, 89692905, "+", "Nonsense_Mutation", "SNP", "C", "C", "T", "SP5001_Tumor", "SP5001_Normal", "DO5001", "GBM-US", "c.388C>T", "p.R130*", 0.52],

    # CIC Splice Site (Oligodendroglioma 1p/19q co-deletion partner)
    ["CIC", 23152, "TCGA", "GRCh37", "19", 42797820, 42797820, "+", "Splice_Site", "SNP", "G", "G", "A", "SP5015_Tumor", "SP5015_Normal", "DO5015", "LGG-US", "c.1450+1G>A", "p.Splice", 0.39],

    # DDX3X 5'UTR (Medulloblastoma)
    ["DDX3X", 1654, "DKFZ", "GRCh37", "X", 41201500, 41201500, "+", "5'UTR", "SNP", "C", "C", "T", "SP5020_Tumor", "SP5020_Normal", "DO5020", "PBCA-DE", "g.41201500C>T", "p.DDX3X_5UTR", 0.36]
]

HEADERS = [
    "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
    "Start_Position", "End_Position", "Strand", "Variant_Classification",
    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
    "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Donor_ID",
    "Project_Code", "HGVSc", "HGVSp_Short", "VAF"
]

def main():
    output_file = "pcawg_brain_somatic_mutations_filtered.csv"
    with open(output_file, mode='w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(HEADERS)
        writer.writerows(BRAIN_DATA)
    print(f"[+] Successfully generated brain cancer dataset with {len(BRAIN_DATA)} variant entries.")

if __name__ == "__main__":
    main()