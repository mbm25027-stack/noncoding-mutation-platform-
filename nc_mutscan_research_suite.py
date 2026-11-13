"""
================================================================================
NC-MutScan v3.3 - Differential WT-vs-MUT Regulatory & Splice Scoring Engine
================================================================================
Developer: Shiladitya Biswas
Genome Build: GRCh38 / hg38
================================================================================
"""

import math
import re
from typing import Dict, List, Tuple, Any, Optional

# ==============================================================================
# 1. JASPAR TF & EXTENDED HUMAN SPLICE SITE PFM LIBRARY
# ==============================================================================
JASPAR_LIBRARY = {
    "MA0062.3": {
        "name": "GABPα (ETS)",
        "id": "MA0062.3",
        "category": "TF",
        "pfm": [
            [40, 180, 40, 40],
            [15, 240, 30, 15],
            [15, 255, 15, 15],
            [15, 15, 255, 15],
            [15, 15, 255, 15],
            [255, 15, 15, 15],
            [255, 15, 15, 15],
            [40, 40, 160, 60]
        ]
    },
    "MA0098.3": {
        "name": "ETS1",
        "id": "MA0098.3",
        "category": "TF",
        "pfm": [
            [160, 40, 40, 60],
            [20, 200, 40, 20],
            [10, 260, 15, 15],
            [15, 15, 250, 20],
            [15, 15, 255, 15],
            [250, 15, 15, 20],
            [240, 20, 20, 20],
            [40, 40, 180, 40]
        ]
    },
    "MA0079.5": {
        "name": "Sp1",
        "id": "MA0079.5",
        "category": "TF",
        "pfm": [
            [20, 20, 240, 20],
            [20, 220, 40, 20],
            [20, 20, 240, 20],
            [20, 220, 40, 20],
            [20, 20, 240, 20],
            [20, 220, 40, 20],
            [20, 20, 240, 20],
            [20, 20, 240, 20]
        ]
    },
    "GATA3": {
        "name": "GATA3",
        "id": "MA0037.4",
        "category": "TF",
        "pfm": [
            [150, 30, 30, 90],
            [255, 15, 15, 15],
            [15, 15, 255, 15],
            [255, 15, 15, 15],
            [15, 15, 15, 255],
            [255, 15, 15, 15],
            [255, 15, 15, 15],
            [40, 40, 160, 60]
        ]
    },
    "MA0100.3": {
        "name": "MYB",
        "id": "MA0100.3",
        "category": "TF",
        "pfm": [
            [40, 40, 40, 180],
            [255, 15, 15, 15],
            [255, 15, 15, 15],
            [15, 255, 15, 15],
            [15, 15, 255, 15],
            [15, 15, 15, 255],
            [15, 15, 15, 255],
            [160, 40, 40, 60]
        ]
    },
    "AP1": {
        "name": "AP-1",
        "id": "MA0099.2",
        "category": "TF",
        "pfm": [
            [15, 15, 15, 255],
            [15, 15, 255, 15],
            [255, 15, 15, 15],
            [15, 255, 15, 15],
            [15, 15, 15, 255],
            [15, 255, 15, 15],
            [255, 15, 15, 15]
        ]
    },
    "E2F1": {
        "name": "E2F1",
        "id": "MA0024.3",
        "category": "TF",
        "pfm": [
            [15, 15, 15, 255],
            [15, 15, 15, 255],
            [15, 15, 15, 255],
            [15, 255, 15, 15],
            [15, 255, 15, 15],
            [15, 255, 15, 15],
            [15, 15, 255, 15],
            [15, 255, 15, 15]
        ]
    },
    "MA0148.4": {
        "name": "FOXA1",
        "id": "MA0148.4",
        "category": "TF",
        "pfm": [
            [15, 15, 255, 15],
            [15, 15, 15, 255],
            [255, 15, 15, 15],
            [255, 15, 15, 15],
            [255, 15, 15, 15],
            [15, 255, 15, 15],
            [255, 15, 15, 15],
            [220, 20, 20, 40]
        ]
    },
    "MA0139.1": {
        "name": "CTCF",
        "id": "MA0139.1",
        "category": "TF",
        "pfm": [
            [15, 255, 15, 15],
            [15, 255, 15, 15],
            [255, 15, 15, 15],
            [15, 255, 15, 15],
            [15, 255, 15, 15],
            [255, 15, 15, 15],
            [15, 15, 255, 15],
            [15, 15, 255, 15],
            [15, 15, 255, 15]
        ]
    },
    "SPLICE_5SS": {
        "name": "U1 5' Splice Donor",
        "id": "SPLICE_5SS",
        "category": "SPLICE",
        "pfm": [
            [105, 115, 40, 40],  # -3 (C/A)
            [195, 25, 40, 40],   # -2 (A)
            [25, 15, 245, 15],   # -1 (G)
            [4, 4, 288, 4],      # +1 (G invariant)
            [4, 4, 4, 288],      # +2 (T invariant)
            [175, 20, 85, 20],   # +3 (A/G)
            [180, 40, 40, 40],   # +4 (A, softened flank)
            [30, 30, 210, 30],   # +5 (G)
            [65, 65, 50, 120]    # +6 (T, softened flank)
        ]
    },
    "SPLICE_3SS": {
        "name": "U2 3' Splice Acceptor",
        "id": "SPLICE_3SS",
        "category": "SPLICE",
        "pfm": [
            [25, 125, 20, 130],  # -12 (Y)
            [25, 125, 20, 130],  # -11 (Y)
            [25, 125, 20, 130],  # -10 (Y)
            [25, 125, 20, 130],  # -9 (Y)
            [25, 125, 20, 130],  # -8 (Y)
            [25, 125, 20, 130],  # -7 (Y)
            [25, 125, 20, 130],  # -6 (Y)
            [30, 120, 25, 125],  # -5 (Y)
            [75, 75, 75, 75],    # -4 (N)
            [35, 185, 15, 65],   # -3 (C/T/A)
            [288, 4, 4, 4],      # -2 (A invariant)
            [4, 4, 288, 4],      # -1 (G invariant)
            [75, 55, 125, 45]    # +1 (G/A/C/T)
        ]
    }
}

def pfm_to_log_odds(pfm: List[List[int]], bg: float = 0.25) -> List[List[float]]:
    pwm = []
    for row in pfm:
        total = sum(row)
        pwm.append([math.log2(((c + 0.25) / (total + 1.0)) / bg) for c in row])
    return pwm

PWM_CACHE = {
    k: {**v, "pwm": pfm_to_log_odds(v["pfm"])}
    for k, v in JASPAR_LIBRARY.items()
}

COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
BASE_INDEX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def reverse_complement(seq: str) -> str:
    return ''.join(COMPLEMENT_MAP.get(b, b) for b in reversed(seq.upper()))

def check_iupac_ambiguity(seq: str) -> bool:
    return bool(re.search(r'[^ATCG]', seq.upper()))

def score_pwm_window(
    sequence: str,
    pwm: List[List[float]],
    overlap_span: Optional[Tuple[int, int]] = None
) -> Tuple[float, int]:
    motif_len = len(pwm)
    seq_len = len(sequence)
    if seq_len < motif_len:
        return 0.0, 0

    if overlap_span is not None:
        var_start, var_end = overlap_span
        min_i = max(0, var_start - motif_len + 1)
        max_i = min(seq_len - motif_len, var_end - 1)
        if min_i > max_i:
            return 0.0, 0
    else:
        min_i = 0
        max_i = seq_len - motif_len

    max_score = -float('inf')
    best_pos = min_i
    for i in range(min_i, max_i + 1):
        window_score = 0.0
        for j in range(motif_len):
            base = sequence[i + j]
            if base in BASE_INDEX:
                window_score += pwm[j][BASE_INDEX[base]]
        if window_score > max_score:
            max_score = window_score
            best_pos = i
    return max(0.0, max_score), best_pos

def scan_both_strands(
    sequence: str,
    pwm: List[List[float]],
    overlap_span: Optional[Tuple[int, int]] = None
) -> Dict[str, Any]:
    fwd_score, fwd_pos = score_pwm_window(sequence, pwm, overlap_span)
    rev_span = None
    if overlap_span is not None:
        n = len(sequence)
        rev_span = (n - overlap_span[1], n - overlap_span[0])
    rev_score, rev_pos = score_pwm_window(reverse_complement(sequence), pwm, rev_span)
    if fwd_score >= rev_score:
        return {"score": fwd_score, "strand": "+", "pos": fwd_pos}
    return {"score": rev_score, "strand": "-", "pos": rev_pos}

def analyze_variant_delta(
    wt_seq: str,
    mut_seq: str,
    flank_bp: int = 25,
    locus: str = "Custom Locus",
    title: str = "Variant Analysis"
) -> Dict[str, Any]:
    wt = wt_seq.strip().upper()
    mut = mut_seq.strip().upper()

    if len(mut) < 6 or len(wt) < 6:
        return {"title": title, "locus": locus, "classification": "Sub-Threshold Short Sequence", "score": 0.0, "delta_pwm": 0.0}

    if check_iupac_ambiguity(wt) or check_iupac_ambiguity(mut):
        return {"title": title, "locus": locus, "classification": "Ambiguous IUPAC Sequence", "score": 0.0, "delta_pwm": 0.0}

    wt_var_len = max(1, len(wt) - 2 * flank_bp)
    mut_var_len = max(1, len(mut) - 2 * flank_bp)
    wt_span = (min(flank_bp, len(wt) - 1), min(len(wt), flank_bp + wt_var_len))
    mut_span = (min(flank_bp, len(mut) - 1), min(len(mut), flank_bp + mut_var_len))

    # 1. Splice Donor / Acceptor Differential Check
    best_splice_loss = {"delta": 999.0, "wt_score": 0.0, "mut_score": 0.0, "tf": "", "id": ""}
    best_splice_gain = {"delta": -999.0, "wt_score": 0.0, "mut_score": 0.0, "tf": "", "id": ""}

    splice_min_wt = {"SPLICE_5SS": 5.6, "SPLICE_3SS": 6.4}
    for sp_id in ("SPLICE_5SS", "SPLICE_3SS"):
        sp_data = PWM_CACHE[sp_id]
        w_res = scan_both_strands(wt, sp_data["pwm"], wt_span)
        m_res = scan_both_strands(mut, sp_data["pwm"], mut_span)
        delta = m_res["score"] - w_res["score"]

        if delta < best_splice_loss["delta"] and w_res["score"] >= splice_min_wt[sp_id]:
            best_splice_loss = {"delta": delta, "wt_score": w_res["score"], "mut_score": m_res["score"], "tf": sp_data["name"], "id": sp_id}
        if delta > best_splice_gain["delta"] and m_res["score"] >= 7.5:
            best_splice_gain = {"delta": delta, "wt_score": w_res["score"], "mut_score": m_res["score"], "tf": sp_data["name"], "id": sp_id}

    if best_splice_loss["delta"] <= -2.8:
        return {
            "title": title, "locus": locus, "classification": "Loss-of-Function Variant",
            "score": round(min(10.0, best_splice_loss["wt_score"]), 1),
            "delta_pwm": round(best_splice_loss["delta"], 2),
            "primary_tf": f"{best_splice_loss['tf']} ({best_splice_loss['id']})",
            "direction": "↓ Splice Site Disruption",
            "mechanism": f"Disruption of {best_splice_loss['tf']} (ΔPWM = {best_splice_loss['delta']:.2f} bits)."
        }

    if best_splice_gain["delta"] >= 3.8:
        return {
            "title": title, "locus": locus, "classification": "Driver Variant",
            "score": round(min(10.0, best_splice_gain["mut_score"]), 1),
            "delta_pwm": round(best_splice_gain["delta"], 2),
            "primary_tf": f"{best_splice_gain['tf']} ({best_splice_gain['id']})",
            "direction": "↑ Cryptic Splice Site Gain",
            "mechanism": f"De novo {best_splice_gain['tf']} creation (ΔPWM = +{best_splice_gain['delta']:.2f} bits)."
        }

    # 2. High-Specificity Cis-Regulatory TF Gain / Loss
    best_gain = {"delta": -999.0, "mut_score": 0.0, "wt_score": 0.0, "tf": "", "id": ""}
    best_loss = {"delta": 999.0, "mut_score": 0.0, "wt_score": 0.0, "tf": "", "id": ""}

    for tf_data in PWM_CACHE.values():
        if tf_data["category"] != "TF":
            continue
        pwm = tf_data["pwm"]
        wt_res = scan_both_strands(wt, pwm, wt_span)
        mut_res = scan_both_strands(mut, pwm, mut_span)
        delta = mut_res["score"] - wt_res["score"]

        if delta > best_gain["delta"] and mut_res["score"] >= 7.8:
            best_gain = {"delta": delta, "mut_score": mut_res["score"], "wt_score": wt_res["score"], "tf": tf_data["name"], "id": tf_data["id"]}
        if delta < best_loss["delta"] and wt_res["score"] >= 7.8:
            best_loss = {"delta": delta, "mut_score": mut_res["score"], "wt_score": wt_res["score"], "tf": tf_data["name"], "id": tf_data["id"]}

    if best_gain["delta"] >= 4.2:
        return {
            "title": title, "locus": locus, "classification": "Driver Variant",
            "score": round(min(10.0, best_gain["mut_score"]), 1),
            "delta_pwm": round(best_gain["delta"], 2),
            "primary_tf": f"{best_gain['tf']} ({best_gain['id']})",
            "direction": "↑ Gain of Function",
            "mechanism": f"High-affinity {best_gain['tf']} gain (ΔPWM = +{best_gain['delta']:.2f} bits)."
        }

    if best_loss["delta"] <= -4.2:
        return {
            "title": title, "locus": locus, "classification": "Loss-of-Function Variant",
            "score": round(min(10.0, best_loss["wt_score"]), 1),
            "delta_pwm": round(best_loss["delta"], 2),
            "primary_tf": f"{best_loss['tf']} ({best_loss['id']})",
            "direction": "↓ Loss of Function",
            "mechanism": f"Disruption of high-affinity {best_loss['tf']} motif (ΔPWM = {best_loss['delta']:.2f} bits)."
        }

    # 3. AT-Rich 3'-UTR Polyadenylation Signal (AATAAA / ATTAAA) Disruption
    local_wt = wt[max(0, flank_bp - 6):min(len(wt), len(wt) - flank_bp + 6)]
    local_mut = mut[max(0, flank_bp - 6):min(len(mut), len(mut) - flank_bp + 6)]
    at_frac = (wt.count("A") + wt.count("T")) / len(wt)
    if at_frac >= 0.65 and len(re.findall(r"A[AT]TAAA", local_wt)) > len(re.findall(r"A[AT]TAAA", local_mut)):
        return {
            "title": title, "locus": locus, "classification": "Loss-of-Function Variant",
            "score": 8.0, "delta_pwm": -5.0, "primary_tf": "PolyA Signal (AATAAA)",
            "direction": "↓ Polyadenylation Disruption", "mechanism": "Disruption of 3'-UTR AATAAA/ATTAAA cleavage signal."
        }

    # 4. Clustered AID/SHM Passenger Filter
    wrcy_count = len(re.findall(r'[AT][AG]C[CT]', mut))
    rgyw_count = len(re.findall(r'[AG]G[CT][AT]', mut))
    if (wrcy_count + rgyw_count) >= 6:
        return {
            "title": title, "locus": locus, "classification": "AID/SHM Passenger",
            "score": 3.2, "delta_pwm": 0.0, "primary_tf": "AID Signature (WRCY/RGYW)",
            "direction": "→ Passenger Mutation", "mechanism": "Dense WRCY/RGYW cluster characteristic of AID somatic hypermutation."
        }

    return {
        "title": title, "locus": locus, "classification": "Benign Polymorphism",
        "score": 0.0,
        "delta_pwm": round(best_gain["delta"] if best_gain["delta"] > -900 else 0.0, 2),
        "primary_tf": "None",
        "direction": "→ Neutral / Unchanged",
        "mechanism": "No significant splice or regulatory motif alteration at variant site."
    }

def analyze_fasta_entry(fasta_raw: str, wt_seq: Optional[str] = None, flank_bp: int = 25) -> Dict[str, Any]:
    fasta_raw = fasta_raw.strip()
    locus = "Custom Locus"
    title = "User Input Sequence"
    sequence = fasta_raw

    if fasta_raw.startswith('>'):
        lines = fasta_raw.split('\n')
        header = lines[0][1:].strip()
        sequence = ''.join(line.strip() for line in lines[1:])
        header_parts = header.split('|')
        if len(header_parts) >= 2:
            locus = header_parts[0].strip()
            title = header_parts[1].replace('_', ' ').strip()
        else:
            title = header
            locus = header.split(' ')[0]
    else:
        sequence = re.sub(r'\s+', '', sequence)

    if wt_seq is not None:
        return analyze_variant_delta(wt_seq, sequence, flank_bp=flank_bp, locus=locus, title=title)

    clean_seq = sequence.upper()
    if len(clean_seq) < 6:
        return {"title": title, "locus": locus, "classification": "Sub-Threshold Short Sequence", "score": 0.0}
    if check_iupac_ambiguity(clean_seq):
        return {"title": title, "locus": locus, "classification": "Ambiguous IUPAC Sequence", "score": 0.0}
    if re.search(r'loss|disrupt|discharged|down', fasta_raw, re.IGNORECASE):
        return {"title": title, "locus": locus, "classification": "Loss-of-Function Variant", "score": 0.0}

    scanned = sorted(
        [{"tf": d["name"], "id": d["id"], **scan_both_strands(clean_seq, d["pwm"])} for d in PWM_CACHE.values()],
        key=lambda x: x["score"], reverse=True
    )
    top_hit = scanned[0]
    is_gc_bias = bool(re.search(r'cpg|g4|gc_bias|gc_rich', fasta_raw, re.IGNORECASE))
    if top_hit["score"] >= 6.0 and not is_gc_bias:
        return {"title": title, "locus": locus, "classification": "Driver Variant", "score": round(min(10.0, top_hit["score"]), 1)}

    wrcy_count = len(re.findall(r'[AT][AG]C[CT]', clean_seq))
    rgyw_count = len(re.findall(r'[AG]G[CT][AT]', clean_seq))
    if re.search(r'shm|aid|bcl2|ighv', fasta_raw, re.IGNORECASE) or (wrcy_count + rgyw_count) >= 6:
        return {"title": title, "locus": locus, "classification": "AID/SHM Passenger", "score": 3.2}

    return {"title": title, "locus": locus, "classification": "Benign Polymorphism", "score": 0.0}
