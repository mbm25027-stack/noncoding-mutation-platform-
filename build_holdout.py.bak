#!/usr/bin/env python3
"""
build_holdout.py - Build a fresh, disjoint, strictly noncoding ClinVar holdout for NC-MutScan.

Reproduces the manuscript's filters (GRCh38, <=15 bp, P/LP vs B/LB, no coding consequence in any
isoform) and adds: review-status filtering, exclusion of every variant used in earlier sets,
many seeded balanced draws, and a manifest with the exact pool counts (the "N = [X]" in Sec. 3.3).

Usage example:
  python build_holdout.py --clinvar clinvar_20261201.vcf.gz --out holdout_v34 \
      --chroms 12-22 --min-stars 2 --n-per-class 75 --seed 42 --n-draws 100 \
      --exclude locked_blind_chr12_22.vcf dev_chr1_11.vcf
"""
import argparse, gzip, hashlib, json, os, random, re, sys
from collections import Counter

CODING_MC = {"missense_variant", "nonsense_variant", "synonymous_variant", "frameshift_variant",
             "inframe_insertion", "inframe_deletion", "stop_lost", "stop_retained_variant",
             "start_lost", "initiator_codon_variant", "coding_sequence_variant",
             "protein_altering_variant", "incomplete_terminal_codon_variant"}
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_single_submitter": 1, "criteria_provided,_conflicting_classifications": 0,
         "criteria_provided,_conflicting_interpretations": 0, "no_assertion_criteria_provided": 0,
         "no_classification_provided": 0, "no_assertion_provided": 0}
PATH = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
BEN = {"Benign", "Likely_benign", "Benign/Likely_benign"}

def opener(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def norm_chrom(c): return c[3:] if c.lower().startswith("chr") else c
def parse_chroms(spec):
    out = set()
    for part in spec.split(","):
        m = re.fullmatch(r"(\d+)-(\d+)", part)
        if m: out.update(str(i) for i in range(int(m[1]), int(m[2]) + 1))
        else: out.add(norm_chrom(part))
    return out
def info_dict(s):
    d = {}
    for kv in s.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v
    return d
def key(c, p, r, a): return f"{norm_chrom(c)}:{p}:{r}:{a}"

def load_exclusions(paths):
    ex = set()
    for p in paths:
        with opener(p) as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                t = line.split("\t")
                ex.add(key(t[0], t[1], t[3], t[4].split(",")[0]))
    return ex

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--clinvar", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--chroms", required=True, help="e.g. 12-22 or 1-11,X")
    ap.add_argument("--max-len", type=int, default=15)
    ap.add_argument("--min-stars", type=int, default=1)
    ap.add_argument("--snv-only", action="store_true")
    ap.add_argument("--n-per-class", type=int, default=75)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-draws", type=int, default=0, help="extra independent balanced draws (seed+i)")
    ap.add_argument("--exclude", nargs="*", default=[], help="VCFs whose variants must not appear")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    chroms = parse_chroms(a.chroms); excl = load_exclusions(a.exclude)
    drops = Counter(); header = []; pools = {"P": [], "B": []}; file_date = None
    with opener(a.clinvar) as f:
        for line in f:
            if line.startswith("##"):
                header.append(line)
                if line.startswith("##fileDate"): file_date = line.strip().split("=", 1)[1]
                continue
            if line.startswith("#"): header.append(line); continue
            t = line.rstrip("\n").split("\t")
            if norm_chrom(t[0]) not in chroms: drops["chrom"] += 1; continue
            ref, alts = t[3], t[4].split(",")
            if len(alts) != 1: drops["multiallelic"] += 1; continue
            alt = alts[0]
            if max(len(ref), len(alt)) > a.max_len: drops["too_long"] += 1; continue
            if a.snv_only and not (len(ref) == 1 and len(alt) == 1): drops["not_snv"] += 1; continue
            info = info_dict(t[7]); sig = info.get("CLNSIG", "")
            if sig in PATH: cls = "P"
            elif sig in BEN: cls = "B"
            else: drops["not_P_or_B"] += 1; continue
            if STARS.get(info.get("CLNREVSTAT", ""), 0) < a.min_stars: drops["review_status"] += 1; continue
            mc = info.get("MC", "")
            if not mc: drops["no_molecular_consequence"] += 1; continue
            terms = {x.split("|")[-1] for x in mc.split(",")}
            if terms & CODING_MC: drops["coding_consequence"] += 1; continue
            if key(t[0], t[1], ref, alt) in excl: drops["previously_used"] += 1; continue
            pools[cls].append(line)
    counts = {"pathogenic_pool": len(pools["P"]), "benign_pool": len(pools["B"])}
    n = min(a.n_per_class, counts["pathogenic_pool"], counts["benign_pool"])
    if n < a.n_per_class: print(f"WARNING: only {n} per class available", file=sys.stderr)

    def write(path, lines):
        with open(path, "w") as o: o.writelines(header + lines)
    write(f"{a.out}/pool_all.vcf", pools["P"] + pools["B"])
    def draw(seed):
        r = random.Random(seed)
        return r.sample(pools["P"], n) + r.sample(pools["B"], n)
    write(f"{a.out}/holdout_balanced.vcf", draw(a.seed))
    for i in range(a.n_draws):
        write(f"{a.out}/draw_{i:03d}.vcf", draw(a.seed + 1 + i))
    man = {"clinvar_file": os.path.basename(a.clinvar), "clinvar_sha256": sha256(a.clinvar),
           "clinvar_fileDate": file_date, "chroms": a.chroms, "min_stars": a.min_stars,
           "max_len": a.max_len, "snv_only": a.snv_only, "seed": a.seed, "n_per_class": n,
           "n_extra_draws": a.n_draws, "excluded_variants_loaded": len(excl),
           "pool_counts_after_filters": counts, "dropped_by_reason": dict(drops),
           "outputs": {fn: sha256(f"{a.out}/{fn}") for fn in sorted(os.listdir(a.out)) if fn.endswith(".vcf")}}
    json.dump(man, open(f"{a.out}/manifest.json", "w"), indent=2)
    print(json.dumps({k: man[k] for k in ("pool_counts_after_filters", "dropped_by_reason", "n_per_class")}, indent=2))
    print(f"Manifest: {a.out}/manifest.json  (chr pool N = {sum(counts.values())})")

if __name__ == "__main__":
    main()