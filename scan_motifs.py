import os
import math
import csv

# ===========================================================================
# 1. DYNAMIC PATH RESOLUTION & AUTOMATIC SEARCH
# ===========================================================================
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def locate_single_file(filename="MA0139.1.pfm"):
    """Locates a specific motif file across potential project directories."""
    candidate_paths = [
        os.path.join(BASE_DIR, "data", "motifs", filename),
        os.path.join(BASE_DIR, "data", filename),
        os.path.join(BASE_DIR, filename)
    ]
    for path in candidate_paths:
        if os.path.isfile(path):
            return path
    return None

def locate_motif_directory():
    """Finds the folder containing JASPAR .jaspar or .pfm files."""
    possible_paths = [
        os.path.join(BASE_DIR, "data", "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar"),
        os.path.join(BASE_DIR, "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar")
    ]
    for path in possible_paths:
        if os.path.isdir(path):
            files = [f for f in os.listdir(path) if f.endswith(('.jaspar', '.pfm'))]
            if files:
                return path
    return None

# Target Locus & Variant Information
TARGET_CHROM = "chr11"
TARGET_START = 5226490
TARGET_END = 5226500

REF_SEQUENCE = "GCTGCCGCCCTGGTGGCTCCACCAGGGGGCAGCTCCCCTGCCCTCGGGGGCTCTGGCCCTGCCCTCCCGGGCCCCCTGGCCCC"
ALT_SEQUENCE = "GCTGCCGCCCTGGTGGCTCCACCAGGGGGCAGTTCCCCTGCCCTCGGGGGCTCTGGCCCTGCCCTCCCGGGCCCCCTGGCCCC"  # C>T Mutation


# ===========================================================================
# 2. ENCODE cCRE BED ANNOTATOR (PURE PYTHON)
# ===========================================================================
def annotate_genomic_ccre(chrom, start_pos, end_pos):
    """
    Scans local ENCODE BED files to label if coordinates overlap with a 
    Promoter-like (PLS), Enhancer-like (ELS), or CTCF element.
    """
    data_dir = os.path.join(BASE_DIR, "data")
    bed_files = [
        ("PLS (Promoter-like)", os.path.join(data_dir, "GRCh38-cCREs.PLS.bed")),
        ("ELS (Enhancer-like)", os.path.join(data_dir, "GRCh38-cCREs.ELS.bed")),
        ("CTCF (Insulator)", os.path.join(data_dir, "GRCh38-cCREs.CA-CTCF.bed")),
        ("cCRE Element", os.path.join(data_dir, "encode_cCREs_GRCh38_combined.bed"))
    ]

    chrom_str = str(chrom) if str(chrom).startswith("chr") else f"chr{chrom}"
    annotations = []

    for label, filepath in bed_files:
        if not os.path.isfile(filepath):
            continue
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#") or line.startswith("track") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    
                    b_chrom, b_start, b_end = parts[0], int(parts[1]), int(parts[2])
                    
                    # Coordinate intersection check
                    if b_chrom == chrom_str and max(start_pos, b_start) < min(end_pos, b_end):
                        accession = parts[3] if len(parts) > 3 else label
                        annotations.append(f"{label} [{accession}]")
                        break  # Found match for this file category
        except Exception:
            continue

    return annotations if annotations else ["Intergenic / Non-cCRE Region"]


# ===========================================================================
# 3. PURE-PYTHON PARSER & LOG-ODDS PWM ENGINE
# ===========================================================================
def reverse_complement(sequence):
    """Computes reverse complement of DNA sequence for antisense (-) strand scanning."""
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(trans)[::-1]

def parse_jaspar_pfm(file_path):
    """Parses raw PFM matrix files (.pfm / .jaspar format)."""
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            lines = [line.strip() for line in f if line.strip()]
    except Exception:
        return None, None, None

    if not lines:
        return None, None, None

    header = lines[0]
    matrix_id = "UNKNOWN"
    tf_name = "UNKNOWN"

    if header.startswith(">"):
        parts = header[1:].split()
        matrix_id = parts[0]
        if len(parts) > 1:
            tf_name = parts[1]
        lines = lines[1:]

    counts = {'A': [], 'C': [], 'G': [], 'T': []}
    bases = ['A', 'C', 'G', 'T']
    
    labeled = any(lines[0].startswith(b) for b in bases)

    try:
        if labeled:
            for line in lines:
                for b in bases:
                    if line.startswith(b):
                        clean_line = line[1:].replace('[', '').replace(']', '').strip()
                        counts[b] = [float(val) for val in clean_line.split()]
        else:
            if len(lines) >= 4:
                for idx, b in enumerate(bases):
                    counts[b] = [float(val) for val in lines[idx].split()]
    except ValueError:
        return None, None, None

    if not counts['A'] or len(counts['A']) == 0:
        return None, None, None

    return matrix_id, tf_name, counts

def build_pwm_log_odds(counts, pseudocount=0.5, background=0.25):
    """Converts raw position count matrix to log-odds Position Weight Matrix."""
    motif_length = len(counts['A'])
    pwm = {'A': [], 'C': [], 'G': [], 'T': []}

    for col in range(motif_length):
        column_total = sum(counts[base][col] for base in ['A', 'C', 'G', 'T'])
        for base in ['A', 'C', 'G', 'T']:
            prob = (counts[base][col] + pseudocount) / (column_total + (4 * pseudocount))
            log_odds = math.log2(prob / background)
            pwm[base].append(log_odds)

    return pwm, motif_length

def max_sliding_pwm_score(pwm, motif_length, target_seq):
    """Calculates peak log-odds PWM score across sequence using sliding window."""
    seq = target_seq.upper()
    if len(seq) < motif_length:
        return -999.0

    max_score = -999.0
    for i in range(len(seq) - motif_length + 1):
        window = seq[i:i + motif_length]
        score = 0.0
        valid_kmer = True
        
        for idx, base in enumerate(window):
            if base in pwm:
                score += pwm[base][idx]
            else:
                valid_kmer = False
                break
                
        if valid_kmer and score > max_score:
            max_score = score

    return max_score

def analyze_variant_impact(file_path):
    """Evaluates sense (+) and antisense (-) strand impact for a single motif file."""
    matrix_id, tf_name, counts = parse_jaspar_pfm(file_path)
    if counts is None:
        return None

    pwm, motif_len = build_pwm_log_odds(counts)

    ref_plus = max_sliding_pwm_score(pwm, motif_len, REF_SEQUENCE)
    alt_plus = max_sliding_pwm_score(pwm, motif_len, ALT_SEQUENCE)
    ref_minus = max_sliding_pwm_score(pwm, motif_len, reverse_complement(REF_SEQUENCE))
    alt_minus = max_sliding_pwm_score(pwm, motif_len, reverse_complement(ALT_SEQUENCE))

    ref_peak = max(ref_plus, ref_minus)
    alt_peak = max(alt_plus, alt_minus)
    best_strand = "+" if alt_plus >= alt_minus else "-"
    delta_pwm = alt_peak - ref_peak

    return {
        "matrix_id": matrix_id,
        "name": tf_name,
        "ref_score": round(ref_peak, 2),
        "alt_score": round(alt_peak, 2),
        "delta_pwm": round(delta_pwm, 2),
        "strand": best_strand
    }


# ===========================================================================
# 4. PIPELINE EXECUTOR WITH ENCODE ANNOTATION
# ===========================================================================
def run_pipeline(delta_threshold=3.0):
    print("=" * 72)
    print("ENCODE cCRE INTEGRATED MOTIF VARIANT SCANNER")
    print("=" * 72)

    # 1. Fetch ENCODE Genomic Context
    ccre_context = annotate_genomic_ccre(TARGET_CHROM, TARGET_START, TARGET_END)
    print(f"[+] Locus           : {TARGET_CHROM}:{TARGET_START}-{TARGET_END}")
    print(f"[+] Genomic Context : {', '.join(ccre_context)}\n")

    # 2. Batch Scan Motifs
    motif_dir = locate_motif_directory()
    if not motif_dir:
        print("[-] Motif dataset directory not found in 'data/'.")
        return

    files = [f for f in os.listdir(motif_dir) if f.endswith(('.jaspar', '.pfm'))]
    print(f"[+] Scanning {len(files)} JASPAR motifs...")

    results = []
    for fname in files:
        fpath = os.path.join(motif_dir, fname)
        res = analyze_variant_impact(fpath)
        if res and abs(res['delta_pwm']) >= delta_threshold:
            res['genomic_context'] = "; ".join(ccre_context)
            results.append(res)

    results.sort(key=lambda x: abs(x['delta_pwm']), reverse=True)

    print(f"\n{'Matrix ID':<12} {'TF Name':<15} {'Strand':<8} {'Ref Score':<10} {'Alt Score':<10} {'Delta PWM':<10}")
    print("-" * 72)
    for r in results[:15]:  # Display top 15 results
        print(f"{r['matrix_id']:<12} {r['name']:<15} {r['strand']:<8} {r['ref_score']:<10} {r['alt_score']:<10} {r['delta_pwm']:<10}")

    # 3. Export to CSV automatically
    csv_file = os.path.join(BASE_DIR, "variant_scan_results.csv")
    keys = ["matrix_id", "name", "strand", "ref_score", "alt_score", "delta_pwm", "genomic_context"]
    with open(csv_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)

    print(f"\n[+] Full batch results with ENCODE context saved to '{os.path.basename(csv_file)}'.")

if __name__ == "__main__":
    run_pipeline(delta_threshold=3.0)