import re
import math
from collections import Counter

KNOWN_DRIVER_MOTIFS = [
    {
        "motif_id": "TERT_PROMOTER_C228T_C250T",
        "pattern": r"TTCCGGGAAC|CCCCTTCCGG",
        "gene": "TERT",
        "score": 8.5,
        "chrom": "chr5",
        "pos": 1295228,
        "cancer_types": ["Melanoma", "Glioblastoma", "Liver / HCC", "Bladder"],
        "encode_ccres": ["PLS"],
        "description": "De novo ETS transcription factor binding motif creation (C228T/C250T promoter mutations)"
    },
    {
        "motif_id": "CTCF_INSULATOR_BOUNDARY_MOTIF",
        "pattern": r"CCACCAGGGGGCGCT",
        "gene": "CTCF_INSULATOR",
        "score": 7.5,
        "chrom": "chr17",
        "pos": 40000000,
        "cancer_types": ["Glioma", "T-ALL", "Gastrointestinal", "Pan-Cancer"],
        "encode_ccres": ["CTCF"],
        "description": "CTCF zinc-finger chromatin boundary / insulator motif disruption"
    },
    {
        "motif_id": "PLEKHS1_PROMOTER_HOTSPOT",
        "pattern": r"GCTAGCCAGGAAGTG",
        "gene": "PLEKHS1",
        "score": 7.0,
        "chrom": "chr10",
        "pos": 115511584,
        "cancer_types": ["Bladder", "Breast", "Head and Neck"],
        "encode_ccres": ["pCRE"],
        "description": "Recurrent non-coding promoter driver hotspot in PLEKHS1"
    },
    {
        "motif_id": "WDR74_PROMOTER_HOTSPOT",
        "pattern": r"CGCGGAAGGCCCCT",
        "gene": "WDR74",
        "score": 7.0,
        "chrom": "chr11",
        "pos": 62604085,
        "cancer_types": ["Melanoma", "Colorectal"],
        "encode_ccres": ["pCRE"],
        "description": "Recurrent non-coding promoter driver motif in WDR74"
    }
]

COMPILED_MOTIFS = [
    {**m, "regex": re.compile(m["pattern"])} for m in KNOWN_DRIVER_MOTIFS
]

def calculate_shannon_entropy(seq: str) -> float:
    if not seq:
        return 0.0
    counts = Counter(seq.upper())
    length = len(seq)
    return -sum((count / length) * math.log2(count / length) for count in counts.values())

def align_and_analyze_sequence(raw_sequence: str, noncoding_db: dict = None) -> dict:
    seq = raw_sequence.strip().upper()
    seq_len = len(seq)
    
    if not seq:
        return {"error": "Empty sequence provided."}
    
    gc_count = seq.count("G") + seq.count("C")
    gc_pct = round((gc_count / seq_len) * 100, 2)
    entropy = calculate_shannon_entropy(seq)
    
    # Filter low-complexity runs and high-GC decoys
    if entropy < 1.2 or (gc_pct > 88.0 and not any(m["regex"].search(seq) for m in COMPILED_MOTIFS)):
        return {
            "status": "Benign / Low-Complexity Repeat",
            "is_driver": False,
            "score": 0.0,
            "sequence_length_bp": seq_len,
            "gc_content_pct": gc_pct,
            "detected_motifs": [],
            "message": "Sequence flagged as low-complexity repeat or non-driver GC decoy."
        }

    detected_hits = []
    max_score = 0.0
    matched_genes = set()
    cancer_types = set()
    encode_ccres = set()

    for item in COMPILED_MOTIFS:
        match = item["regex"].search(seq)
        if match:
            max_score = max(max_score, item["score"])
            matched_genes.add(item["gene"])
            cancer_types.update(item["cancer_types"])
            encode_ccres.update(item["encode_ccres"])
            detected_hits.append({
                "motif_id": item["motif_id"],
                "matched_subsequence": match.group(0),
                "start_in_seq": match.start(),
                "end_in_seq": match.end(),
                "gene": item["gene"],
                "chrom": item["chrom"],
                "pos": item["pos"],
                "cancer_types": item["cancer_types"],
                "encode_ccres": item["encode_ccres"],
                "description": item["description"]
            })

    is_driver = len(detected_hits) > 0 and max_score >= 6.0

    return {
        "status": "Pathogenic / Non-Coding Driver Sequence" if is_driver else "Benign / Safe Sequence",
        "is_driver": is_driver,
        "score": max_score if is_driver else 0.0,
        "gene": ", ".join(matched_genes) if matched_genes else "None",
        "sequence_length_bp": seq_len,
        "gc_content_pct": gc_pct,
        "detected_motifs": detected_hits,
        "cancer_types": sorted(list(cancer_types)),
        "encode_ccres": sorted(list(encode_ccres)),
        "message": f"Identified {len(detected_hits)} pathogenic non-coding driver motif(s)." if is_driver else "No driver motifs detected."
    }