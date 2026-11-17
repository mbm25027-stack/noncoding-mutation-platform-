#!/usr/bin/env python3
"""
extract_spliceai_scores.py - stream SpliceAI scores out of a LARGE source and keep only
the variants in your pool. Standard library only. Nothing is unzipped onto disk.

Sources it understands
  * a .zip holding MANY per-gene tables (.tsv/.tsv.gz), e.g. AbSplice_DNA_hg38_snvs.zip
  * a single table (.tsv/.csv, optionally .gz or a one-member .zip)
  * an Illumina SpliceAI VCF (--vcf), plain or .gz

Output: TSV  key<TAB>ds_max   (key = chrom:pos:ref:alt, chrom without "chr").
If a variant appears several times (several genes/tissues/rows) the maximum is kept.

STEP 1 - look (prints the zip layout, then the header and first rows of the first gene file):
  python extract_spliceai_scores.py --peek --scores-file AbSplice_DNA_hg38_snvs.zip

STEP 2 - extract (columns are auto-detected from the header: chrom/pos/ref/alt and a column
named delta_score; override with --chrom-col --pos-col --ref-col --alt-col --score-col).
Recommended: --gtf (same GENCODE GTF you used for the junction analysis) so that only the gene
files that overlap your pool variants are read - much faster:
  python extract_spliceai_scores.py --scores-file AbSplice_DNA_hg38_snvs.zip ^
     --pool check_blind\\pool_all.vcf --gtf $gtf --out spliceai_scores.tsv
Use the splice-inclusive pool VCF (check_blind\\pool_all.vcf) so one file serves both analyses.
Variants that never appear in the output are "unscorable" in the analysis.
"""
import argparse, bisect, gzip, io, os, sys, time, zipfile
from collections import defaultdict

DATA_EXT = (".tsv", ".tsv.gz", ".csv", ".csv.gz", ".txt", ".txt.gz")
ALIASES = dict(chrom=("chrom", "chr", "chromosome", "#chrom", "chromosome_name"),
               pos=("pos", "position", "start", "pos_hg38"),
               ref=("ref", "reference", "ref_allele"),
               alt=("alt", "alternate", "alt_allele"))


def wrap(raw, name):
    if name.lower().endswith(".gz"):
        raw = gzip.open(raw)
    return io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")


def data_members(z):
    names = [n for n in z.namelist() if not n.endswith("/")]
    sel = [n for n in names if n.lower().endswith(DATA_EXT)]
    return sel or names


def open_plain(path):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def pool_keys(pool_path):
    keys = set()
    with open(pool_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split("\t", 5)
            keys.add(f"{t[0].replace('chr', '')}:{t[1]}:{t[3].upper()}:{t[4].upper()}")
    return keys


def genes_overlapping_pool(gtf_path, keys):
    pos = defaultdict(list)
    for k in keys:
        c, p, _ = k.split(":", 2)
        pos[c].append(int(p))
    for c in pos:
        pos[c].sort()
    out = set()
    with open_plain(gtf_path) as g:
        for line in g:
            if line[0] == "#":
                continue
            t = line.split("\t", 9)
            if len(t) < 9 or t[2] != "gene":
                continue
            c = t[0][3:] if t[0].startswith("chr") else t[0]
            L = pos.get(c)
            if not L:
                continue
            i = bisect.bisect_left(L, int(t[3]))
            if i < len(L) and L[i] <= int(t[4]):
                a = t[8]
                j = a.find('gene_id "')
                if j >= 0:
                    out.add(a[j + 9:a.find('"', j + 9)].split(".")[0])
    return out


def sniff_delim(line):
    return "\t" if line.count("\t") >= line.count(",") else ","


def find_col(header, spec, alias_key):
    if spec:
        if spec.isdigit():
            return int(spec)
        if spec not in header:
            raise SystemExit(f"Column '{spec}' not found. Header columns: {header}")
        return header.index(spec)
    low = [h.lower() for h in header]
    for a in ALIASES[alias_key]:
        if a in low:
            return low.index(a)
    raise SystemExit(f"Could not find a '{alias_key}' column automatically. Header columns: {header}\n"
                     f"Pass it explicitly, e.g. --{alias_key}-col NAME")


def find_score_cols(header, spec):
    if spec:
        return [int(x) if x.strip().isdigit() else header.index(x.strip()) for x in spec.split(",")]
    low = [h.lower() for h in header]
    if "delta_score" in low:
        return [low.index("delta_score")]
    cand = [i for i, h in enumerate(low) if "delta" in h and "score" in h]
    if cand:
        return cand
    raise SystemExit(f"Could not find the score column automatically. Header columns: {header}\n"
                     f"Pass it explicitly with --score-col NAME")


def scan_table(f, keys, best, a, cache):
    """Scan one delimited table (own header). Returns number of lines read."""
    first = next(f, "")
    if not first:
        return 0
    delim = sniff_delim(first)
    header = [h.strip().strip('"') for h in first.rstrip("\r\n").split(delim)]
    ck = tuple(header)
    if ck not in cache:
        cache[ck] = (find_col(header, a.chrom_col, "chrom"), find_col(header, a.pos_col, "pos"),
                     find_col(header, a.ref_col, "ref"), find_col(header, a.alt_col, "alt"),
                     find_score_cols(header, a.score_col))
    ic, ip, ir, ia, isc = cache[ck]
    need = max([ic, ip, ir, ia] + isc)
    n = 0
    for line in f:
        n += 1
        t = line.rstrip("\r\n").split(delim)
        if len(t) <= need:
            continue
        chrom = t[ic].strip('"')
        chrom = chrom[3:] if chrom.startswith("chr") else chrom
        k = f"{chrom}:{t[ip].strip(chr(34))}:{t[ir].strip(chr(34)).upper()}:{t[ia].strip(chr(34)).upper()}"
        if k not in keys:
            continue
        try:
            v = max(float(t[j].strip('"')) for j in isc if t[j] not in ("", "NA", "NaN", "nan"))
        except ValueError:
            continue
        if k not in best or v > best[k]:
            best[k] = v
    return n


def scan_vcf(f, keys, best, t0):
    n = 0
    for line in f:
        n += 1
        if line[0] == "#":
            continue
        t = line.split("\t", 8)
        if len(t) < 8:
            continue
        chrom = t[0][3:] if t[0].startswith("chr") else t[0]
        k0 = f"{chrom}:{t[1]}:{t[3].upper()}:"
        for alt in t[4].upper().split(","):
            k = k0 + alt
            if k not in keys:
                continue
            info = t[7]
            i = info.find("SpliceAI=")
            if i < 0:
                continue
            end = info.find(";", i)
            for ent in info[i + 9:(end if end >= 0 else None)].split(","):
                p = ent.split("|")
                if len(p) < 6 or p[0].upper() != alt:
                    continue
                try:
                    v = max(float(x) for x in p[2:6] if x not in (".", ""))
                except ValueError:
                    continue
                if k not in best or v > best[k]:
                    best[k] = v
        if n % 20_000_000 == 0:
            print(f"  {n / 1e6:,.0f}M lines, {len(best):,} found ({time.time() - t0:,.0f}s)")
    return n


def write_out(best, path):
    tmp = path + ".tmp"
    with open(tmp, "w") as o:
        o.write("key\tds_max\n")
        for k, v in best.items():
            o.write(f"{k}\t{v}\n")
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores-file", required=True)
    ap.add_argument("--member", default=None, help="one member name inside a .zip (default: all data members)")
    ap.add_argument("--peek", action="store_true")
    ap.add_argument("--vcf", action="store_true", help="Illumina SpliceAI VCF (INFO SpliceAI=...)")
    ap.add_argument("--pool", help="pool VCF whose variants to keep")
    ap.add_argument("--gtf", help="GENCODE GTF: read only gene files overlapping pool variants (zip of per-gene files)")
    ap.add_argument("--out", default="spliceai_scores.tsv")
    ap.add_argument("--chrom-col"); ap.add_argument("--pos-col"); ap.add_argument("--ref-col")
    ap.add_argument("--alt-col"); ap.add_argument("--score-col", help="column name(s), comma separated -> max")
    a = ap.parse_args()

    is_zip = a.scores_file.lower().endswith(".zip")
    z = zipfile.ZipFile(a.scores_file) if is_zip else None
    members = data_members(z) if is_zip else []
    if a.member:
        members = [a.member]

    if a.peek:
        if is_zip:
            tot = sum(z.getinfo(m).compress_size for m in members) / 1e9
            print(f"Zip has {len(members):,} data members, {tot:.1f} GB compressed. First members:")
            for m in members[:5]:
                i = z.getinfo(m)
                print(f"  {m}  ({i.compress_size / 1e6:.1f} MB compressed, {i.file_size / 1e6:.1f} MB uncompressed)")
            f = wrap(z.open(members[0]), members[0])
            print(f"\nHeader and first rows of {members[0]}:")
        else:
            f = open_plain(a.scores_file)
        for i, line in enumerate(f):
            print(line.rstrip("\r\n")[:700])
            if i >= 5:
                break
        return

    if not a.pool:
        raise SystemExit("--pool is required (the pool VCF whose variants to keep)")
    keys = pool_keys(a.pool)
    print(f"{len(keys):,} pool variants to look up")
    best, t0, nlines, cache = {}, time.time(), 0, {}

    if is_zip and len(members) > 1 and not a.vcf:
        if a.gtf:
            want = genes_overlapping_pool(a.gtf, keys)
            have = {os.path.basename(m).split(".")[0]: m for m in members}
            present = sorted(g for g in want if g in have)
            print(f"GTF genes overlapping pool variants: {len(want):,}; present as files in the zip: {len(present):,}; "
                  f"absent from the zip: {len(want) - len(present):,} (variants only in absent genes will count as unscorable)")
            members = [have[g] for g in present]
        print(f"Reading {len(members):,} gene files ({sum(z.getinfo(m).compress_size for m in members) / 1e9:.1f} GB compressed)")
        for i, m in enumerate(members, 1):
            nlines += scan_table(wrap(z.open(m), m), keys, best, a, cache)
            if i % 100 == 0 or i == len(members):
                el = time.time() - t0
                print(f"  {i:,}/{len(members):,} files | {nlines / 1e6:,.1f}M rows | {len(best):,} pool variants found | "
                      f"{el / 60:,.1f} min elapsed, ~{el / i * (len(members) - i) / 60:,.0f} min left")
            if i % 500 == 0:
                write_out(best, a.out)
    else:
        if is_zip:
            m = members[0]
            f = wrap(z.open(m), m)
        else:
            f = open_plain(a.scores_file)
        nlines = scan_vcf(f, keys, best, t0) if a.vcf else scan_table(f, keys, best, a, cache)

    write_out(best, a.out)
    print(f"Read {nlines:,} rows in {(time.time() - t0) / 60:,.1f} min; kept {len(best):,} of {len(keys):,} pool variants "
          f"({len(best) / max(1, len(keys)):.1%}) -> {a.out}")


if __name__ == "__main__":
    main()
