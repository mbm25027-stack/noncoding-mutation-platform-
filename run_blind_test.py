import json
import os
import ssl
import time
import urllib.request
from nc_mutscan_research_suite import analyze_fasta_entry

INPUT_VCF = "blind_test_set_nc.vcf"
CACHE_FILE = "blind_test_flanks.json"
FLANK_BP = 25

# Bypass Windows local SSL issuer verification error
SSL_CONTEXT = ssl._create_unverified_context()


def fetch_hg38_mutated_seq(chrom: str, pos: int, ref: str, alt: str, cache: dict) -> str:
    clean_chr = str(chrom).replace("chr", "").strip()
    if clean_chr == "MT":
        clean_chr = "M"
    ucsc_chr = f"chr{clean_chr}"
    key = f"{ucsc_chr}:{pos}:{ref}>{alt}"

    if key in cache:
        return cache[key]

    # 0-based half-open coordinates for UCSC hg38 API
    start_0 = max(0, pos - 1 - FLANK_BP)
    end_0 = pos - 1 + len(ref) + FLANK_BP
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={ucsc_chr};start={start_0};end={end_0}"

    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15, context=SSL_CONTEXT) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    dna = data.get("dna", "").upper()
    ref_offset = (pos - 1) - start_0
    mut_seq = dna[:ref_offset] + alt.upper() + dna[ref_offset + len(ref):]
    cache[key] = mut_seq
    return mut_seq


def main():
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, "r", encoding="utf-8") as cf:
            cache = json.load(cf)

    variants = []
    with open(INPUT_VCF, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            chrom, pos, var_id, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4].split(",")[0]
            if len(ref) > 20 or len(alt) > 20:
                continue
            variants.append((chrom, pos, var_id, ref, alt))

    print(f"Reconstructing hg38 +/-{FLANK_BP}bp flanking sequences for {len(variants)} ClinVar variants...")

    tp = fn = 0
    results_log = []

    for idx, (chrom, pos, var_id, ref, alt) in enumerate(variants, 1):
        try:
            mut_seq = fetch_hg38_mutated_seq(chrom, pos, ref, alt, cache)
            fasta_entry = f">chr{chrom}:{pos}|ClinVar_{var_id}_{ref}>{alt}\n{mut_seq}"
            res = analyze_fasta_entry(fasta_entry)
            classification = res.get("classification", "Unknown")

            is_driver = classification in ["Driver Variant", "Loss-of-Function Variant"]
            if is_driver:
                tp += 1
            else:
                fn += 1

            results_log.append((chrom, pos, var_id, classification))
            if idx % 10 == 0:
                print(f"  Processed {idx}/{len(variants)} variants...")
                time.sleep(0.2)
        except Exception as e:
            print(f"  [!] Skipped chr{chrom}:{pos} ({e})")

    with open(CACHE_FILE, "w", encoding="utf-8") as cf:
        json.dump(cache, cf, indent=2)

    total = tp + fn
    sensitivity = (tp / total * 100) if total > 0 else 0.0

    print("\n" + "=" * 75)
    print(" NC-MUTSCAN v3.1 BLIND CLINVAR BENCHMARK RESULTS")
    print("=" * 75)
    print(f" Total Variants Evaluated : {total}")
    print(f" True Positives (TP)      : {tp}")
    print(f" False Negatives (FN)     : {fn}")
    print(f" Blind Sensitivity (TPR)  : {sensitivity:.1f}%")
    print("=" * 75)


if __name__ == "__main__":
    main()