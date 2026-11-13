import json
import os
import ssl
import time
import urllib.request
from nc_mutscan_research_suite import analyze_fasta_entry

DEV_VCF = "dev_chr1_11.vcf"
CACHE_FILE = "dev_flanks.json"
FLANK_BP = 25
SSL_CONTEXT = ssl._create_unverified_context()


def fetch_wt_and_mut_seq(chrom: str, pos: int, ref: str, alt: str, cache: dict):
    clean_chr = str(chrom).replace("chr", "").strip()
    ucsc_chr = f"chr{clean_chr}"
    key = f"{ucsc_chr}:{pos}:{ref}>{alt}"

    if key in cache:
        return cache[key]["wt"], cache[key]["mut"]

    start_0 = max(0, pos - 1 - FLANK_BP)
    end_0 = pos - 1 + len(ref) + FLANK_BP
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={ucsc_chr};start={start_0};end={end_0}"

    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15, context=SSL_CONTEXT) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    dna = data.get("dna", "").upper()
    ref_offset = (pos - 1) - start_0
    wt_seq = dna
    mut_seq = dna[:ref_offset] + alt.upper() + dna[ref_offset + len(ref):]

    cache[key] = {"wt": wt_seq, "mut": mut_seq}
    return wt_seq, mut_seq


def main():
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, "r", encoding="utf-8") as cf:
            cache = json.load(cf)

    records = []
    with open(DEV_VCF, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            label = parts[0]
            chrom, pos, var_id, ref, alt = parts[1], int(parts[2]), parts[3], parts[4], parts[5].split(",")[0]
            records.append((label, chrom, pos, var_id, ref, alt))

    tp = fn = tn = fp = 0
    false_negatives = []
    false_positives = []

    print(f"Auditing {len(records)} dev variants (chr1-chr11)...")
    for idx, (label, chrom, pos, var_id, ref, alt) in enumerate(records, 1):
        try:
            wt_seq, mut_seq = fetch_wt_and_mut_seq(chrom, pos, ref, alt, cache)
            fasta_entry = f">chr{chrom}:{pos}|{label}_{var_id}_{ref}>{alt}\n{mut_seq}"
            res = analyze_fasta_entry(fasta_entry)
            pred_class = res.get("classification", "Unknown")
            is_driver = pred_class in ["Driver Variant", "Loss-of-Function Variant"]

            if label == "PATHOGENIC":
                if is_driver:
                    tp += 1
                else:
                    fn += 1
                    false_negatives.append((f"chr{chrom}:{pos} ({ref}>{alt})", pred_class, res.get("score", 0), mut_seq))
            else:
                if is_driver:
                    fp += 1
                    false_positives.append((f"chr{chrom}:{pos} ({ref}>{alt})", pred_class, res.get("score", 0), mut_seq))
                else:
                    tn += 1

            if idx % 15 == 0:
                print(f"  Processed {idx}/{len(records)}...")
                time.sleep(0.15)
        except Exception as e:
            print(f"  [!] Skipped chr{chrom}:{pos}: {e}")

    with open(CACHE_FILE, "w", encoding="utf-8") as cf:
        json.dump(cache, cf, indent=2)

    sens = (tp / (tp + fn) * 100) if (tp + fn) else 0.0
    spec = (tn / (tn + fp) * 100) if (tn + fp) else 0.0

    print("\n" + "=" * 75)
    print(" DEV SET (CHR1-CHR11) BASELINE AUDIT")
    print("=" * 75)
    print(f" TP: {tp} | FN: {fn} | TN: {tn} | FP: {fp}")
    print(f" Dev Sensitivity (TPR): {sens:.1f}%")
    print(f" Dev Specificity (TNR): {spec:.1f}%")
    print("-" * 75)
    print(" SAMPLE FALSE NEGATIVES (Missed Pathogenic):")
    for item in false_negatives[:8]:
        print(f"   {item[0]} -> {item[1]} | score={item[2]} | seq={item[3]}")
    print(" SAMPLE FALSE POSITIVES (Benign Flagged as Driver):")
    for item in false_positives[:8]:
        print(f"   {item[0]} -> {item[1]} | score={item[2]} | seq={item[3]}")
    print("=" * 75)


if __name__ == "__main__":
    main()