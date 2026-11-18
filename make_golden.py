#!/usr/bin/env python3
"""
make_golden.py -- run ONCE, locally, on the verified frozen build, to record the expected outputs.
The golden file holds the flank sequences, so CI needs no genome FASTA.

  python tests/make_golden.py --fasta "<GRCh38.d1.vd1.fa>" --dev dev_chr1_11.vcf \
         --blind locked_blind_chr12_22_clean --blind-n 50 --out tests/golden/golden_v3.3.0.json

Needs spliceVarDB_eval.py (Fasta, build_flanks, load_engine) in the project folder.
"""
import argparse, json, os, random, re, sys
sys.path.insert(0, os.getcwd())
from spliceVarDB_eval import Fasta, resolve_fasta, build_flanks, load_engine, norm_chrom, open_text, FLANK_BP

NEG = ("benign", "neutral", "polymorphism", "passenger", "sub-threshold", "subthreshold", "non-driver", "compet")


def parse_variants(path):
    """Standard VCF, or the label-first layout (LABEL, chrom, pos, id, ref, alt, ...). Returns dicts."""
    out = []
    with open_text(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            for off in (0, 1):
                if len(t) < off + 5:
                    continue
                try:
                    chrom, pos = norm_chrom(t[off]), int(t[off + 1])
                except ValueError:
                    continue
                ref, alt = t[off + 3].strip().upper(), t[off + 4].strip().upper()
                if not (re.fullmatch(r"[ACGTN]+", ref) and re.fullmatch(r"[ACGTN]+", alt)):
                    continue
                text = (t[0] if off == 1 else (t[7] if len(t) > 7 else "")).lower()
                if "clnsig=" in text:
                    text = text.split("clnsig=")[1].split(";")[0]
                label = 0 if any(n in text for n in ("benign",)) else (1 if "pathogenic" in text else None)
                out.append({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "label": label})
                break
    return out


def pick(rows, n, seed):
    rng = random.Random(seed)
    pos = [r for r in rows if r["label"] == 1]
    neg = [r for r in rows if r["label"] == 0]
    if pos and neg:
        h = n // 2
        return rng.sample(pos, min(h, len(pos))) + rng.sample(neg, min(n - h, len(neg)))
    return rng.sample(rows, min(n, len(rows)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--dev", required=True)
    ap.add_argument("--blind", required=True)
    ap.add_argument("--blind-n", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20261118)
    ap.add_argument("--release", default="RELEASE.json")
    ap.add_argument("--server-dir", default=".")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rel = json.load(open(a.release, encoding="utf-8"))
    score, server = load_engine(a.server_dir)
    short, exp = getattr(server, "CORE_SHORT_SHA256", None), getattr(server, "EXPECTED_BUILD_PREFIX", None)
    if not short or short != exp or not rel["core_sha256"].startswith(short):
        raise SystemExit(f"Build not verified (engine {short}, expected {exp}, pinned {rel['core_sha256'][:12]}). Stop.")

    fa = Fasta(resolve_fasta(a.fasta))
    sets = [("dev", parse_variants(a.dev)), ("blind", pick(parse_variants(a.blind), a.blind_n, a.seed))]
    records, skipped = [], 0
    for name, rows in sets:
        for r in rows:
            fl, err = build_flanks(fa, r["chrom"], r["pos"], r["ref"], r["alt"])
            if err:
                skipped += 1
                continue
            s, call, res = score(*fl)
            if res is None:
                skipped += 1
                continue
            records.append({"id": f"{name}:{r['chrom']}:{r['pos']}:{r['ref']}>{r['alt']}", "set": name,
                            "label": r["label"], "wt_seq": fl[0], "mut_seq": fl[1],
                            "S_raw": s, "call": call, "classification": str(res.get("classification", ""))})
    meta = {"release": rel["git_tag"], "core_sha256": rel["core_sha256"], "flank_bp": FLANK_BP,
            "generated_by": "make_golden.py", "blind_seed": a.seed, "blind_n": a.blind_n,
            "skipped": skipped, "counts": {k: sum(r["set"] == k for r in records) for k in ("dev", "blind")}}
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump({"meta": meta, "records": records}, f, indent=1)
    print(f"Wrote {a.out}: {meta['counts']}, skipped {skipped}")


if __name__ == "__main__":
    main()
