import scan_motifs

def align_and_analyze_sequence(raw_sequence: str, noncoding_db: dict = None) -> dict:
    """
    Main sequence analysis function used by server.py.
    Evaluates sequence entropy, GC percentage, and driver motifs.
    """
    seq = raw_sequence.strip().upper()
    seq_len = len(seq)
    
    if not seq:
        return {"error": "Empty sequence provided."}
    
    gc_count = seq.count("G") + seq.count("C")
    gc_pct = round((gc_count / seq_len) * 100, 2)
    entropy = scan_motifs.calculate_shannon_entropy(seq)
    
    # 1. Low-Complexity & Mononucleotide Repeat Suppression
    if entropy < 1.25:
        return {
            "status": "Benign / Low-Complexity Repeat",
            "is_driver": False,
            "score": 0.0,
            "sequence_length_bp": seq_len,
            "gc_content_pct": gc_pct,
            "shannon_entropy": entropy,
            "detected_motifs": [],
            "message": "Sequence flagged as low-complexity repeat (filtered by entropy)."
        }
    
    # 2. Motif Scanning
    detected_hits = scan_motifs.scan_sequence_motifs(seq)
    
    # 3. High-GC Decoy Filtering (rejects GC-dense sequences lacking explicit driver motifs)
    if gc_pct > 80.0 and not detected_hits:
        return {
            "status": "Benign / High-GC Background Decoy",
            "is_driver": False,
            "score": 0.0,
            "sequence_length_bp": seq_len,
            "gc_content_pct": gc_pct,
            "shannon_entropy": entropy,
            "detected_motifs": [],
            "message": "High-GC promoter background sequence without validated driver consensus."
        }
    
    max_score = max([hit["score"] for hit in detected_hits], default=0.0)
    matched_genes = list({hit["gene"] for hit in detected_hits})
    cancer_types = sorted(list({ct for hit in detected_hits for ct in hit["cancer_types"]}))
    encode_ccres = sorted(list({ec for hit in detected_hits for ec in hit["encode_ccres"]}))
    
    is_driver = len(detected_hits) > 0 and max_score >= 6.0
    
    return {
        "status": "Pathogenic / Non-Coding Driver Sequence" if is_driver else "Benign / Safe Sequence",
        "is_driver": is_driver,
        "score": max_score if is_driver else 0.0,
        "gene": ", ".join(matched_genes) if matched_genes else "None",
        "sequence_length_bp": seq_len,
        "gc_content_pct": gc_pct,
        "shannon_entropy": entropy,
        "detected_motifs": detected_hits,
        "cancer_types": cancer_types,
        "encode_ccres": encode_ccres,
        "message": f"Identified {len(detected_hits)} pathogenic non-coding driver motif(s)." if is_driver else "No driver motifs detected."
    }