"""
run_item3_exact.py
Re-evaluates the 15-case verification set, 50-variant held-out set, and
54-locus GRCh38 benchmark using NC-MutScan v3.3 (analyze_fasta_entry).
"""

import os
import re
from nc_mutscan_research_suite import analyze_fasta_entry
import nc_mutscan_heldout_benchmark as hb

POSITIVE_CLASSES = {"Driver Variant", "Loss-of-Function Variant"}

# 1. 50-Variant Held-Out Set under v3.3 (single-sequence mode)
tp = fp = tn = fn = 0
for item in hb.HELD_OUT_DATASET:
    res = analyze_fasta_entry(item["fasta"], wt_seq=None)
    pred_pos = res.get("classification") in POSITIVE_CLASSES
    exp_pos = item["label"] == "Driver"
    if exp_pos and pred_pos:
        tp += 1
    elif exp_pos and not pred_pos:
        fn += 1
    elif not exp_pos and pred_pos:
        fp += 1
    else:
        tn += 1

print("=== ITEM 3A: 50-Variant Held-Out Benchmark under NC-MutScan v3.3 ===")
print(f"  N = {tp+fp+tn+fn} | TP={tp}, FP={fp}, TN={tn}, FN={fn}")
print(f"  Sensitivity = {tp/(tp+fn)*100:.2f}% | Specificity = {tn/(tn+fp)*100:.2f}% | "
      f"Precision = {tp/max(1,tp+fp)*100:.2f}% | Accuracy = {(tp+tn)/(tp+fp+tn+fn)*100:.2f}%")

# 2. 54-Locus GRCh38 Benchmark under v3.3 (single-sequence mode)
BENCHMARK_DRIVERS = [
    (">chr5:1295113|TERT_promoter|-124C>T", "CCCACGTGGCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCCCTTCCT"),
    (">chr10:113824420|PLEKHS1_promoter|C>T", "GGAGGCCGCGGCCGCCGGAACTCGGGCCCGCCGCCCAGCCCCGCCCCCGCC"),
    (">chr11:111957575|SDHD_promoter|C>T", "CCCGCCCCCGCCGCCGGAACTCGGCCGCCGCCCCCGCCGCCCGCCCCCGCC"),
    (">chr1:47210120|TAL1_enhancer|MYB_insert", "AGGCTCGTGGAACCGTTAGCTAGCTAGCGTTAGCTAGGCTAGCTAGCTA"),
]

records_54 = list(BENCHMARK_DRIVERS)
if os.path.exists("benchmark_200.fasta"):
    with open("benchmark_200.fasta", "rt", encoding="utf-8") as f:
        for block in f.read().strip().split(">"):
            if not block.strip():
                continue
            lines = block.strip().splitlines()
            records_54.append((">" + lines[0], "".join(lines[1:]).upper().replace(" ", "")))

counts_54 = {}
driver_calls_54 = []
for hdr, seq in records_54:
    res = analyze_fasta_entry(f"{hdr}\n{seq}", wt_seq=None)
    cls = res.get("classification", "Unknown")
    counts_54[cls] = counts_54.get(cls, 0) + 1
    if cls in POSITIVE_CLASSES:
        driver_calls_54.append((hdr, cls, res.get("top_motif", res.get("motif", ""))))

print(f"\n=== ITEM 3B: {len(records_54)}-Locus GRCh38 Benchmark under NC-MutScan v3.3 ===")
for cls, cnt in sorted(counts_54.items(), key=lambda x: -x[1]):
    print(f"  {cls:28s}: {cnt:2d} / {len(records_54)} ({cnt/len(records_54)*100:.1f}%)")
print("  Driver-called loci:", [x[0] for x in driver_calls_54])

# 3. Check 15-case verification in verify_full_system.py if present
if os.path.exists("verify_full_system.py"):
    txt = open("verify_full_system.py", "rt", encoding="utf-8").read()
    pairs = re.findall(r'\(\s*["\'](>[^"\']+)["\']\s*,\s*["\']([ACGTNacgtn]+)["\']\s*,\s*["\']([^"\']+)["\']\s*\)', txt)
    if pairs:
        matched = 0
        for hdr, seq, exp_cls in pairs:
            res = analyze_fasta_entry(f"{hdr}\n{seq}", wt_seq=None)
            got = res.get("classification", "")
            if (exp_cls in POSITIVE_CLASSES and got in POSITIVE_CLASSES) or (exp_cls not in POSITIVE_CLASSES and got not in POSITIVE_CLASSES):
                matched += 1
        print(f"\n=== ITEM 3C: Verification Set ({len(pairs)} cases) under v3.3 ===")
        print(f"  Concordance: {matched}/{len(pairs)} ({matched/len(pairs)*100:.1f}%)")
