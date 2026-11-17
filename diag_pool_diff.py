#!/usr/bin/env python3
"""Compare the old clean pool (check_blind\\pool_all.vcf minus the evaluation's
excluded consequences) with the rebuilt pool (check_blind_v2\\pool_all.vcf).
Shows what the extra / missing variants look like. Standard library only.

  python diag_pool_diff.py check_blind\\pool_all.vcf check_blind_v2\\pool_all.vcf
"""
import sys
from collections import Counter

OLD_EXCL = {"nonsense", "splice_donor_variant", "splice_acceptor_variant",
            "inframe_indel", "no_sequence_alteration"}


def load(path, excl):
    d = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
            terms = frozenset(x.split("|")[-1] for x in info.get("MC", "").split(",") if x)
            if terms & excl:
                continue
            sig = info.get("CLNSIG", "").lower()
            lab = "P" if ("pathogenic" in sig and "benign" not in sig) else "B"
            key = f"{t[0].replace('chr', '')}:{t[1]}:{t[3].upper()}:{t[4].upper()}"
            d[key] = (lab, terms, len(t[3]), len(t[4]), info.get("CLNSIG", ""),
                      info.get("CLNREVSTAT", ""), info.get("CLNVC", ""))
    return d


def bucket(n):
    return "1 (SNV-size)" if n == 1 else "2-5" if n <= 5 else "6-10" if n <= 10 else "11-15" if n <= 15 else ">15"


def describe(title, d, keys):
    print(f"\n=== {title}: {len(keys):,} variants ===")
    if not keys:
        return
    print("label:", dict(Counter(d[k][0] for k in keys)))
    print("variant type (CLNVC):", dict(Counter(d[k][6] for k in keys).most_common(8)))
    print("max(len ref, len alt):", dict(Counter(bucket(max(d[k][2], d[k][3])) for k in keys)))
    print("consequence terms:", dict(Counter(x for k in keys for x in d[k][1]).most_common(12)))
    print("review status:", dict(Counter(d[k][5] for k in keys).most_common(5)))
    print("CLNSIG:", dict(Counter(d[k][4] for k in keys).most_common(5)))
    print("examples:", list(keys)[:5])


old = load(sys.argv[1], OLD_EXCL)
new = load(sys.argv[2], set())
print(f"old clean pool: {len(old):,} (P {sum(v[0]=='P' for v in old.values()):,})")
print(f"rebuilt pool:   {len(new):,} (P {sum(v[0]=='P' for v in new.values()):,})")
only_new = [k for k in new if k not in old]
only_old = [k for k in old if k not in new]
describe("ONLY in rebuilt pool", new, only_new)
describe("ONLY in old pool", old, only_old)