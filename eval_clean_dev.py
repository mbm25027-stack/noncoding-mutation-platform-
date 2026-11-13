import json
import os
import ssl
import time
import urllib.request
from collections import Counter
from nc_mutscan_research_suite import analyze_variant_delta

DEV_VCF = "dev_chr1_11.vcf"
CACHE_FILE = "dev_flanks_clean.json"
FLANK_BP = 25
SSL_CTX = ssl._create_unverified_context()

def fetch_seqs(chrom: str, pos: int, ref: str, alt: str, cache: dict):
    clean_chr = str(chrom).replace("chr", "").strip()
    ucsc_chr = f"chr{clean_chr}"
    key = f"{ucsc_chr}:{pos}:{ref}>{alt}"
    if key in cache:
        return cache[key]["wt"], cache[key]["mut"]

    start_0 = max(0, pos - 1 - FLANK_BP)
    end_0 = pos - 1 + len(ref) + FLANK_BP
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={ucsc_chr};start={start_0};end={end_0}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15, context=SSL_CTX) as resp:
        dna = json.loads(resp.read().decode("utf-8")).get("dna", "").upper()

    offset = (pos - 1) - start_0
    mut_dna = dna[:offset] + alt.upper() + dna[offset + len(ref):]
    cache[key] = {"wt": dna, "mut": mut_dna}
    return dna, mut_dna

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
        records.append((parts[0], parts[1], int(parts[2]), parts[4], parts[5].split(",")[0], parts[8] if len(parts) > 8 else ""))

print(f"Evaluating {len(records)} strictly non-coding dev variants (chr1-chr11)...")
tp = fn = tn = fp = 0
mc_stats = Counter()
fn_samples, fp_samples = [], []

for idx, (label, chrom, pos, ref, alt, info) in enumerate(records, 1):
    try:
        wt, mut = fetch_seqs(chrom, pos, ref, alt, cache)
        res = analyze_variant_delta(wt, mut, flank_bp=FLANK_BP, locus=f"chr{chrom}:{pos}")
        is_driver = res["classification"] in ["Driver Variant", "Loss-of-Function Variant"]

        mc_tag = "other"
        for field in info.split(";"):
            if field.startswith("MC="):
                mc_tag = field.split("|")[1] if "|" in field else field

        if label == "PATHOGENIC":
            if is_driver:
                tp += 1
                mc_stats[f"TP:{mc_tag}"] += 1
            else:
                fn += 1
                mc_stats[f"FN:{mc_tag}"] += 1
                if len(fn_samples) < 6:
                    fn_samples.append((f"chr{chrom}:{pos} {ref}>{alt}", mc_tag, res.get("delta_pwm", 0), mut))
        else:
            if is_driver:
                fp += 1
                if len(fp_samples) < 6:
                    fp_samples.append((f"chr{chrom}:{pos} {ref}>{alt}", res.get("primary_tf", ""), res.get("delta_pwm", 0)))
            else:
                tn += 1

        if idx % 25 == 0:
            with open(CACHE_FILE, "w", encoding="utf-8") as cf:
                json.dump(cache, cf)
            print(f"  Processed {idx}/{len(records)}...")
            time.sleep(0.1)
    except Exception as e:
        print(f"  [!] Skipped chr{chrom}:{pos}: {e}")

with open(CACHE_FILE, "w", encoding="utf-8") as cf:
    json.dump(cache, cf)

sens = tp / (tp + fn) * 100 if (tp + fn) else 0.0
spec = tn / (tn + fp) * 100 if (tn + fp) else 0.0
prec = tp / (tp + fp) * 100 if (tp + fp) else 0.0
f1 = (2 * tp / (2 * tp + fp + fn)) if (2 * tp + fp + fn) else 0.0

print("\n" + "=" * 75)
print(" CLEAN NON-CODING DEV SET (CHR1-CHR11) v3.2 RESULTS")
print("=" * 75)
print(f" TP: {tp} | FN: {fn} | TN: {tn} | FP: {fp}")
print(f" Dev Sensitivity (TPR) : {sens:.1f}%")
print(f" Dev Specificity (TNR) : {spec:.1f}%")
print(f" Dev Precision (PPV)   : {prec:.1f}%")
print(f" Dev F1-Score          : {f1:.4f}")
print("-" * 75)
print(" CONSEQUENCE BREAKDOWN:", dict(mc_stats))
print(" SAMPLE FNs:", fn_samples)
print(" SAMPLE FPs:", fp_samples)
print("=" * 75)
