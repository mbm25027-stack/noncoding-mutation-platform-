import json
from nc_mutscan_research_suite import analyze_variant_delta

DEV_VCF = "dev_chr1_11.vcf"
CACHE_FILE = "dev_flanks.json"

with open(CACHE_FILE, "r", encoding="utf-8") as cf:
    cache = json.load(cf)

tp = fn = tn = fp = 0
with open(DEV_VCF, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.strip().split("\t")
        label = parts[0]
        chrom, pos, var_id, ref, alt = parts[1], int(parts[2]), parts[3], parts[4], parts[5].split(",")[0]
        key = f"chr{chrom.replace('chr', '').strip()}:{pos}:{ref}>{alt}"
        if key not in cache:
            continue

        res = analyze_variant_delta(cache[key]["wt"], cache[key]["mut"], flank_bp=25, locus=f"chr{chrom}:{pos}")
        is_driver = res["classification"] in ["Driver Variant", "Loss-of-Function Variant"]

        if label == "PATHOGENIC":
            if is_driver:
                tp += 1
            else:
                fn += 1
        else:
            if is_driver:
                fp += 1
            else:
                tn += 1

sens = (tp / (tp + fn) * 100) if (tp + fn) else 0.0
spec = (tn / (tn + fp) * 100) if (tn + fp) else 0.0
prec = (tp / (tp + fp) * 100) if (tp + fp) else 0.0
f1 = (2 * tp / (2 * tp + fp + fn)) if (2 * tp + fp + fn) else 0.0

print("=" * 75)
print(" NC-MUTSCAN v3.2 DEV SET (CHR1-CHR11) EVALUATION")
print("=" * 75)
print(f" TP: {tp} | FN: {fn} | TN: {tn} | FP: {fp}")
print(f" Dev Sensitivity (TPR) : {sens:.1f}%")
print(f" Dev Specificity (TNR) : {spec:.1f}%")
print(f" Dev Precision (PPV)   : {prec:.1f}%")
print(f" Dev F1-Score          : {f1:.4f}")
print("=" * 75)
