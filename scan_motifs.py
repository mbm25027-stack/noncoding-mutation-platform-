import re
import math
from collections import Counter

# Driver Motifs: Pre-defined consensus patterns for non-coding driver hotspots
KNOWN_DRIVER_MOTIFS = [
    {
        "motif_id": "TERT_PROMOTER_C228T_C250T",
        "pattern": r"TTCCGGGAAC|CCCCTTCCGG|CCTTCCGG|TTCCGGGA",
        "gene": "TERT",
        "score": 9.0,
        "chrom": "chr5",
        "pos": 1295228,
        "cancer_types": ["Melanoma", "Glioblastoma", "Liver / HCC", "Bladder"],
        "encode_ccres": ["PLS"],
        "description": "De novo ETS transcription factor binding motif creation (C228T/C250T promoter mutations)"
    },
    {
        "motif_id": "CTCF_INSULATOR_BOUNDARY_MOTIF",
        "pattern": r"CCACCAGGGGGCGCT|CCACCAGGGGG|CAGGGGGCGCT",
        "gene": "CTCF_INSULATOR",
        "score": 8.0,
        "chrom": "chr17",
        "pos": 40000000,
        "cancer_types": ["Glioma", "T-ALL", "Gastrointestinal", "Pan-Cancer"],
        "encode_ccres": ["CTCF"],
        "description": "CTCF zinc-finger chromatin boundary / insulator motif disruption"
    },
    {
        "motif_id": "PLEKHS1_PROMOTER_HOTSPOT",
        "pattern": r"GCTAGCCAGGAAGTG|CCAGGAAGTG",
        "gene": "PLEKHS1",
        "score": 7.5,
        "chrom": "chr10",
        "pos": 115511584,
        "cancer_types": ["Bladder", "Breast", "Head and Neck"],
        "encode_ccres": ["pCRE"],
        "description": "Recurrent non-coding promoter driver hotspot in PLEKHS1"
    },
    {
        "motif_id": "WDR74_PROMOTER_HOTSPOT",
        "pattern": r"CGCGGAAGGCCCCT|CGCGGAAGG",
        "gene": "WDR74",
        "score": 7.5,
        "chrom": "chr11",
        "pos": 62604085,
        "cancer_types": ["Melanoma", "Colorectal"],
        "encode_ccres": ["pCRE"],
        "description": "Recurrent non-coding promoter driver motif in WDR74"
    },
    {
        "motif_id": "SDHD_PROMOTER_HOTSPOT",
        "pattern": r"CGAGCCGCCCCT|GCCCCCGAGCC",
        "gene": "SDHD",
        "score": 7.0,
        "chrom": "chr11",
        "pos": 111959625,
        "cancer_types": ["Melanoma", "Paraganglioma"],
        "encode_ccres": ["pCRE"],
        "description": "Recurrent ETS promoter binding driver mutation in SDHD"
    },
    {
        "motif_id": "FOXA1_ENHANCER_HOTSPOT",
        "pattern": r"GTAAACAAACA|TGTTTGTTTA",
        "gene": "FOXA1",
        "score": 7.5,
        "chrom": "chr14",
        "pos": 38062000,
        "cancer_types": ["Prostate", "Breast"],
        "encode_ccres": ["dCRE"],
        "description": "Recurrent non-coding enhancer driver mutation in FOXA1 regulatory region"
    }
]

# Pre-compile regexes at module load time for < 1ms match throughput
COMPILED_MOTIFS = [
    {**m, "regex": re.compile(m["pattern"], re.IGNORECASE)} for m in KNOWN_DRIVER_MOTIFS
]

def calculate_shannon_entropy(seq: str) -> float:
    """Calculates Shannon entropy of sequence to filter low-complexity repeats."""
    if not seq:
        return 0.0
    seq_upper = seq.upper()
    length = len(seq_upper)
    counts = Counter(seq_upper)
    return round(-sum((c / length) * math.log2(c / length) for c in counts.values()), 4)

def scan_sequence_motifs(seq: str):
    """Scans sequence against compiled driver motifs."""
    seq_clean = seq.strip().upper()
    hits = []
    for item in COMPILED_MOTIFS:
        match = item["regex"].search(seq_clean)
        if match:
            hits.append({
                "motif_id": item["motif_id"],
                "matched_subsequence": match.group(0),
                "start_in_seq": match.start(),
                "end_in_seq": match.end(),
                "score": item["score"],
                "gene": item["gene"],
                "chrom": item["chrom"],
                "pos": item["pos"],
                "cancer_types": item["cancer_types"],
                "encode_ccres": item["encode_ccres"],
                "description": item["description"]
            })
    return hits