import os
from flask import Flask, request, jsonify
from scan_motifs import (
    annotate_genomic_ccre,
    locate_motif_directory,
    analyze_variant_impact,
    BASE_DIR
)

app = Flask(__name__)

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint to verify backend status."""
    return jsonify({
        "status": "online",
        "pipeline": "Strand-Aware Noncoding Mutation Analyzer",
        "jaspar_version": "2024 CORE Vertebrates"
    })

@app.route('/api/analyze', methods=['POST'])
def analyze_variant():
    """
    POST payload format:
    {
        "chrom": "chr11",
        "start": 5226490,
        "end": 5226500,
        "ref_seq": "GCTGCCGCCCTGGTGGCTCCACCAGGGGGCAGCTCCCCTGCCCTCGGGGGCTCTGGCCCTGCCCTCCCGGGCCCCCTGGCCCC",
        "alt_seq": "GCTGCCGCCCTGGTGGCTCCACCAGGGGGCAGTTCCCCTGCCCTCGGGGGCTCTGGCCCTGCCCTCCCGGGCCCCCTGGCCCC",
        "threshold": 3.0
    }
    """
    data = request.json or {}
    
    chrom = data.get("chrom", "chr11")
    start = int(data.get("start", 5226490))
    end = int(data.get("end", 5226500))
    ref_seq = data.get("ref_seq", "")
    alt_seq = data.get("alt_seq", "")
    threshold = float(data.get("threshold", 3.0))

    if not ref_seq or not alt_seq:
        return jsonify({"error": "Missing reference or alternate sequence."}), 400

    # 1. Annotate ENCODE Genomic Context
    ccre_context = annotate_genomic_ccre(chrom, start, end)

    # 2. Run Motif Scan over Local Library
    motif_dir = locate_motif_directory()
    if not motif_dir:
        return jsonify({"error": "JASPAR motif library directory not found."}), 500

    results = []
    files = [f for f in os.listdir(motif_dir) if f.endswith(('.jaspar', '.pfm'))]

    for fname in files:
        fpath = os.path.join(motif_dir, fname)
        res = analyze_variant_impact(fpath)
        if res and abs(res['delta_pwm']) >= threshold:
            results.append(res)

    # Sort by absolute Delta PWM magnitude
    results.sort(key=lambda x: abs(x['delta_pwm']), reverse=True)

    return jsonify({
        "status": "success",
        "locus": f"{chrom}:{start}-{end}",
        "genomic_context": ccre_context,
        "total_motifs_scanned": len(files),
        "significant_hits_count": len(results),
        "results": results
    })

if __name__ == "__main__":
    # Standard Flask dev server on port 5000
    app.run(host="127.0.0.1", port=5000, debug=True)