import json
import ssl
import time
import urllib.request
from nc_mutscan_research_suite import analyze_variant_delta

BLIND_VCF = "locked_blind_chr12_22.vcf"
FLANK_BP = 25
SSL_CTX = ssl._create_unverified_context()

def fetch_seqs(chrom: str, pos: int, ref: str, alt: str):
    clean_chr = str(chrom).replace("chr", "").strip()
    ucsc_chr = f"chr{clean_chr}"
    start_0 = max(0, pos - 1 - FLANK_BP)
    end_0 = pos - 1 + len(ref) + FLANK_BP
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={ucsc_chr};start={start_0};end={end_0}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15, context=SSL_CTX) as resp:
        dna = json.loads(resp.read().decode("utf-8")).get("dna", "").upper()
    offset = (pos - 1) - start_0
    mut_dna = dna[:offset] + alt.upper() + dna[offset + len(ref):]
    return dna, mut_dna

records = []
with open(BLIND_VCF, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        records.append((parts[0], parts[1], int(parts[2]), parts[4], parts[5].split(",")[0]))

print(f"Running ONE-SHOT final evaluation on {len(records)} untouched blind variants (chr12-chr22)...")
tp = fn = tn = fp = 0

for idx, (label, chrom, pos, ref, alt) in enumerate(records, 1):
    try:
        wt, mut = fetch_seqs(chrom, pos, ref, alt)
        res = analyze_variant_delta(wt, mut, flank_bp=FLANK_BP, locus=f"chr{chrom}:{pos}")
        is_driver = res["classification"] in ["Driver Variant", "Loss-of-Function Variant"]

        if label == "PATHOGENIC":
            if is_driver: tp += 1
            else: fn += 1
        else:
            if is_driver: fp += 1
            else: tn += 1

        if idx % 25 == 0:
            print(f"  Evaluated {idx}/{len(records)} blind variants...")
            time.sleep(0.1)
    except Exception as e:
        print(f"  [!] Skipped chr{chrom}:{pos}: {e}")

total = tp + fn + tn + fp
acc = (tp + tn) / total * 100 if total else 0.0
sens = tp / (tp + fn) * 100 if (tp + fn) else 0.0
spec = tn / (tn + fp) * 100 if (tn + fp) else 0.0
prec = tp / (tp + fp) * 100 if (tp + fp) else 0.0
f1 = (2 * tp / (2 * tp + fp + fn)) if (2 * tp + fp + fn) else 0.0

print("\n" + "=" * 80)
print(" NC-MUTSCAN v3.3 LOCKED BLIND BENCHMARK (CHR12-CHR22 | ZERO LEAKAGE)")
print("=" * 80)
print(f" Total Blind Variants Evaluated : {total} (75 Pathogenic / 75 Benign)")
print(f" Confusion Matrix               : TP={tp} | FP={fp} | TN={tn} | FN={fn}")
print(f" Blind Accuracy                 : {acc:.2f}%")
print(f" Blind Sensitivity (TPR)        : {sens:.2f}%")
print(f" Blind Specificity (TNR)        : {spec:.2f}%")
print(f" Blind Precision (PPV)          : {prec:.2f}%")
print(f" Blind F1-Score                 : {f1:.4f}")
print("=" * 80)
