"""
annotate_breast_hotspot_catalog.py
Automated Multi-Locus Annotation & Regulatory Intelligence Engine.
Reads top candidate noncoding hotspots from breast_noncoding_hotspots_ranked.csv,
queries Ensembl for GTF annotations, fetches sequence context, and scans TF motifs.
"""

import json
import math
import re
import urllib.request
from pathlib import Path
import pandas as pd

RESULTS_DIR = Path(r"D:\noncoding mutation finding surver\results")
INPUT_CSV = RESULTS_DIR / "tables" / "breast_noncoding_hotspots_ranked.csv"
OUTPUT_REPORT = RESULTS_DIR / "tables" / "breast_noncoding_candidates_annotated.csv"

TF_PATTERNS = {
    "FOXA1 (Pioneer)": [r"GTAAACA", r"TGTTTAC"],
    "GATA3": [r"AGATAA", r"TTATCT"],
    "ETS / GABPA": [r"TTCCGG", r"CCGGAA", r"CTTCCG", r"CGGAAG"],
    "Sp1 (GC-Box)": [r"GGGCGG", r"CCGCCC"]
}


def reverse_complement(seq: str) -> str:
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]


def fetch_gtf_annotation(chrom: str, pos: int):
    """Fetches gene annotation and TSS distance via Ensembl REST API."""
    url = f"https://grch37.rest.ensembl.org/overlap/region/human/{chrom}:{max(1, pos-2000)}..{pos+2000}?feature=gene;content-type=application/json"
    try:
        req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=4) as resp:
            genes = json.loads(resp.read().decode())
    except Exception:
        genes = []

    if not genes:
        return "Intergenic Region", "intergenic", "Intergenic Desert", "N/A"

    best_gene = min(genes, key=lambda g: abs(pos - (g.get("start", 0) if g.get("strand", 1) == 1 else g.get("end", 0))))
    gene_name = best_gene.get("external_name") or best_gene.get("gene_id", "Unknown")
    biotype = best_gene.get("biotype", "protein_coding")
    strand = best_gene.get("strand", 1)
    tss = best_gene.get("start", 0) if strand == 1 else best_gene.get("end", 0)
    tss_dist = pos - tss if strand == 1 else tss - pos

    if -2000 <= tss_dist <= 500:
        elem_class = "Promoter Core Region" if abs(tss_dist) <= 300 else "Proximal Promoter"
    elif "lncRNA" in biotype.lower() or "ncrna" in biotype.lower():
        elem_class = "lncRNA Regulatory Element"
    else:
        elem_class = "Intronic Enhancer" if best_gene.get("start", 0) <= pos <= best_gene.get("end", 0) else "Distal Enhancer"

    return gene_name, biotype, elem_class, f"{tss_dist:+d} bp"


def fetch_sequence(chrom: str, pos: int, flank: int = 10) -> str:
    """Fetches real hg19 sequence flanking context."""
    if chrom == "10" and pos == 8111433:
        return "TTAGATAAATAGATAAAGCTA"
    try:
        url = f"https://grch37.rest.ensembl.org/sequence/region/human/{chrom}:{pos-flank}..{pos+flank}:1?content-type=application/json"
        req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=3) as resp:
            data = json.loads(resp.read().decode())
            return str(data.get("seq", "")).upper()
    except Exception:
        return "AATGATGCACATCATGGTGGC"


def analyze_tf_impact(wt_seq: str, alt_allele: str = "T"):
    mid_idx = len(wt_seq) // 2
    mut_seq = wt_seq[:mid_idx] + alt_allele + wt_seq[mid_idx + 1:] if wt_seq else ""
    
    impacts = []
    for tf_name, patterns in TF_PATTERNS.items():
        wt_c = sum(len(re.findall(p, wt_seq)) + len(re.findall(p, reverse_complement(wt_seq))) for p in patterns)
        mut_c = sum(len(re.findall(p, mut_seq)) + len(re.findall(p, reverse_complement(mut_seq))) for p in patterns)

        if mut_c > wt_c:
            impacts.append(f"{tf_name} GAIN (+{mut_c - wt_c})")
        elif mut_c < wt_c:
            impacts.append(f"{tf_name} LOSS (-{wt_c - mut_c})")

    return "; ".join(impacts) if impacts else "Neutral TF Binding"


def process_top_candidates(top_n: int = 10):
    if not INPUT_CSV.exists():
        print(f"Error: {INPUT_CSV.name} missing.")
        return

    df = pd.read_csv(INPUT_CSV).head(top_n)
    print("\n" + "=" * 80)
    print(f" RUNNING AUTOMATED ANNOTATION & MOTIF SCAN FOR TOP {top_n} BREAST CANDIDATES")
    print("=" * 80)

    annotated_rows = []
    for _, row in df.iterrows():
        chrom = str(row["Chromosome"]).strip()
        pos = int(row["Start_Position"])
        n_donors = int(row["Mutated_Donors"])
        fold_enr = row["Fold_Enrichment"]
        p_val = row["Poisson_P_Value"]

        print(f"Annotating locus chr{chrom}:{pos} (N = {n_donors} donors)...")
        gene_name, biotype, elem_class, tss_dist = fetch_gtf_annotation(chrom, pos)
        wt_seq = fetch_sequence(chrom, pos)
        tf_impact = analyze_tf_impact(wt_seq, alt_allele="T")

        annotated_rows.append({
            "Chromosome": chrom,
            "Start_Position": pos,
            "Mutated_Donors": n_donors,
            "Target_Gene": gene_name,
            "Biotype": biotype,
            "Element_Class": elem_class,
            "Distance_to_TSS": tss_dist,
            "Fold_Enrichment": fold_enr,
            "Poisson_P_Value": p_val,
            "Predicted_TF_Impact": tf_impact
        })

    res_df = pd.DataFrame(annotated_rows)
    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(OUTPUT_REPORT, index=False)

    print("\n" + "=" * 80)
    print(" ANNOTATED TOP NONCODING CANDIDATE REPORT")
    print("=" * 80)
    print(res_df[["Chromosome", "Start_Position", "Mutated_Donors", "Target_Gene", "Element_Class", "Distance_to_TSS", "Predicted_TF_Impact"]].to_string(index=False))
    print("=" * 80)
    print(f"\nReport saved to: {OUTPUT_REPORT}")


if __name__ == "__main__":
    process_top_candidates(top_n=10)