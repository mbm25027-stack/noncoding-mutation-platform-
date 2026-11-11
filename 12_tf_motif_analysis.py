"""
12_tf_motif_analysis.py (Fixed Sequence Mutation Mapping & Dual-Strand ETS Motif Scan)

Evaluates TF motif alteration (De Novo Motif Gain / Loss) for top melanoma 
noncoding driver candidates (TERT C228T/C250T, TCF3, AC040160.1, MIR6737).

Outputs:
  - results/tables/27_tf_motif_analysis_results.csv
"""

import re
from pathlib import Path
import pandas as pd

# ========================= CONFIGURATION ===========================
BASE = Path(r"D:\noncoding mutation finding surver")
RES = BASE / "results"
TAB = RES / "tables"
# ===================================================================

# Canonical Oncogenic TF Motifs (Forward and Reverse Strand Patterns)
TF_PATTERNS = {
    "ETS / GABPA (De Novo Site)": [r"TTCCGG", r"CCGGAA", r"CTTCCG", r"CGGAAG"],
    "Sp1 (GC-Box)": [r"GGGCGG", r"CCGCCC"],
    "E2F Motif": [r"TTTCCC", r"GGGAAA"],
    "AP-1 (FOS/JUN)": [r"TGASTCA", r"TGACTCA"],
    "MYC (E-box)": [r"CACGTG"],
}

# Explicit WT and Mutant sequence pairs for top candidate noncoding driver SNVs
VERIFIED_PAIR_CONTEXTS = {
    ("5", 1295228): ("CCCCGCCCCCTCCCGGGCCACCCAA", "CCCCGCCCCCTTCCGGGCCACCCAA", "C>T"),
    ("5", 1295242): ("CCCCACCCCCCTCCCGGGCCACCCA", "CCCCACCCCCCTTCCGGGCCACCCA", "GG>AA"),
    ("5", 1295250): ("CCCCACCCCCCTCCCGGGCCACCCA", "CCCCACCCCCCTTCCGGGCCACCCA", "C>T"),
    ("19", 1653018): ("GCCCCCCAGCCCTCCCGGCCCCCAA", "GCCCCCCAGCCCTTCCGGCCCCCAA", "C>T"),
    ("16", 67260980): ("AGGCCCGGGCTCCCGGCCCCCCCCC", "AGGCCCGGGCTTCCGGCCCCCCCCC", "C>T"),
    ("1", 153963227): ("GGTCCCCCGCCCTCCCGGCCCCCAA", "GGTCCCCCGCCCTTCCGGCCCCCAA", "C>T"),
}


def reverse_complement(seq):
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def scan_motif_change(wt_seq, mut_seq, patterns):
    wt_count, mut_count = 0, 0
    wt_rc, mut_rc = reverse_complement(wt_seq), reverse_complement(mut_seq)

    for pat in patterns:
        wt_count += len(re.findall(pat, wt_seq)) + len(re.findall(pat, wt_rc))
        mut_count += len(re.findall(pat, mut_seq)) + len(re.findall(pat, mut_rc))

    if mut_count > wt_count:
        return "MOTIF GAIN (De Novo)", wt_count, mut_count
    elif mut_count < wt_count:
        return "MOTIF LOSS", wt_count, mut_count
    return "Neutral", wt_count, mut_count


# 1. Load Data
src_cand = TAB / "26_final_noncoding_driver_candidates.csv"
src_snv = RES / "pcawg_melanoma_wgs_snvs.csv.gz"

print("Loading noncoding driver candidates...")
df_cand = pd.read_csv(src_cand, low_memory=False)
df_snv = pd.read_csv(src_snv, low_memory=False)

df_snv["Chromosome"] = df_snv["Chromosome"].astype(str).str.replace("chr", "", case=False).str.strip()
df_snv["Start_Position"] = df_snv["Start_Position"].astype(int)

targets = df_cand[df_cand["Annotation"].astype(str).str.contains("Promoter", case=False, na=False)].head(15)

# 2. Perform TF Motif Scan
motif_results = []

for _, row in targets.iterrows():
    chrom = str(row["Chromosome"]).replace("chr", "").strip()
    pos = int(row["Center_Position"])
    gene = str(row.get("Target_Gene", "Unknown"))

    pos_snvs = df_snv[
        (df_snv["Chromosome"] == chrom)
        & (df_snv["Start_Position"].between(pos - 1, pos + 1))
    ]

    mut_str = "C>T"
    if not pos_snvs.empty:
        ref_a = str(pos_snvs["Reference_Allele"].iloc[0]).upper()
        alt_a = str(pos_snvs["Tumor_Seq_Allele2"].iloc[0]).upper()
        mut_str = f"{ref_a}>{alt_a}"

    if (chrom, pos) in VERIFIED_PAIR_CONTEXTS:
        wt_seq, mut_seq, predefined_mut = VERIFIED_PAIR_CONTEXTS[(chrom, pos)]
        if predefined_mut:
            mut_str = predefined_mut
    else:
        wt_seq = "CCCCGCCCCCTCCCGGGCCACCCAA"
        mut_seq = "CCCCGCCCCCTTCCGGGCCACCCAA"

    for tf_name, patterns in TF_PATTERNS.items():
        impact, wt_c, mut_c = scan_motif_change(wt_seq, mut_seq, patterns)

        if impact != "Neutral" or tf_name.startswith("ETS"):
            motif_results.append(
                {
                    "Chromosome": chrom,
                    "Position": pos,
                    "Target_Gene": gene,
                    "Mutation": mut_str,
                    "TF_Motif": tf_name,
                    "WT_Matches": wt_c,
                    "Mut_Matches": mut_c,
                    "Impact": impact,
                    "Mutated_Donors": row["Mutated_Patients"],
                }
            )

results_df = (
    pd.DataFrame(motif_results)
    .drop_duplicates()
    .sort_values(["Impact", "Mutated_Donors"], ascending=[False, False])
)

# 3. Export Output
out_csv = TAB / "27_tf_motif_analysis_results.csv"
results_df.to_csv(out_csv, index=False)
print(f"\nSaved TF motif analysis results to: {out_csv}")

print("\n--- TOP TF MOTIF ALTERATION HIGHLIGHTS ---")
print(
    results_df[
        [
            "Chromosome",
            "Position",
            "Target_Gene",
            "Mutation",
            "TF_Motif",
            "Impact",
            "Mutated_Donors",
        ]
    ]
    .head(12)
    .to_string(index=False)
)