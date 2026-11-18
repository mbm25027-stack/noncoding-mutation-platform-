#!/usr/bin/env python3
"""
spliceVarDB_eval.py -- external test of the FROZEN NC-MutScan v3.3 engine on SpliceVarDB.
No parameter is changed. Needs numpy only.

It copies the protocol of heldout_full_pool_eval.py:
  * score in-process with server.run_v33_prediction({"wt_seq","mut_seq","flank_bp":25})
  * plus-strand GRCh38 window: 25 bp + ref allele + 25 bp, ref must match, no strand flip
  * call = is_frozen_driver_call(classification text); AUROC uses S_raw
  * max allele length 15 (same as the ClinVar pool manifest)
It then removes every variant that is also in the ClinVar calibration set
(dev_chr1_11.vcf) and reports how many that was.

CORRECTION in this version: read_keys() now handles the calibration file that has an
extra leading label column (PATHOGENIC<TAB>chrom<TAB>pos<TAB>id<TAB>ref<TAB>alt ...)
as well as a standard VCF. The old version read the label as the chromosome, so the
overlap removal silently did nothing.

Run from the project folder (the one with server.py):
  python spliceVarDB_eval.py --svdb splicevardb_hg38.tsv --fasta "<path to GRCh38.d1.vd1.fa>" --smoke 200
  python spliceVarDB_eval.py --svdb splicevardb_hg38.tsv --fasta "<path to GRCh38.d1.vd1.fa>" --out eval_splicevardb

If the script can't find a column, it stops and prints the header it saw.
Use --col-chrom / --col-pos / --col-ref / --col-alt / --col-label to name them.
"""
import argparse, csv, gzip, json, math, os, re, sys, time
from collections import Counter
import numpy as np

FLANK_BP = 25


# ------------------------------------------------------------ FASTA (same reader as heldout_full_pool_eval.py)
class Fasta:
    def __init__(self, path):
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai):
            print(f"Building index {fai} (one-time, may take a few minutes)...")
            self._build(fai)
        self.idx = {}
        with open(fai) as f:
            for line in f:
                n, L, off, lb, lw = line.rstrip("\n").split("\t")[:5]
                self.idx[n] = (int(L), int(off), int(lb), int(lw))
        self.fh = open(path, "rb")

    def _build(self, fai_path):
        ents, name, offset = [], None, 0
        length = seq_off = lb = lw = 0
        with open(self.path, "rb") as f:
            for line in f:
                if line.startswith(b">"):
                    if name:
                        ents.append((name, length, seq_off, lb, lw))
                    name = line[1:].split()[0].decode()
                    length, seq_off, lb, lw = 0, offset + len(line), None, None
                else:
                    s = line.rstrip(b"\r\n")
                    if lb is None:
                        lb, lw = len(s), len(line)
                    length += len(s)
                offset += len(line)
            if name:
                ents.append((name, length, seq_off, lb, lw))
        with open(fai_path, "w") as o:
            for e in ents:
                o.write("\t".join(map(str, e)) + "\n")

    def contig(self, chrom):
        for n in (f"chr{chrom}", chrom):
            if n in self.idx:
                return n
        return None

    def fetch(self, name, start, end):  # 0-based, half-open
        L, off, lb, lw = self.idx[name]
        if start < 0 or end > L:
            return None
        first, last = start // lb, (end - 1) // lb
        b0 = off + first * lw + start % lb
        b1 = off + last * lw + (end - 1) % lb + 1
        self.fh.seek(b0)
        return self.fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def resolve_fasta(p):
    if os.path.isfile(p):
        return p
    if os.path.isdir(p):
        c = [os.path.join(r, f) for r, _, fs in os.walk(p) for f in fs
             if f.lower().endswith((".fa", ".fasta"))]
        if c:
            best = max(c, key=os.path.getsize)
            print(f"{p} is a folder; using {best}")
            return best
    raise SystemExit(f"Cannot find a FASTA at {p}")


# ------------------------------------------------------------ frozen engine (copied from the frozen scripts)
def is_frozen_driver_call(res):   # verbatim from build_dev_flanks_and_loco.py
    cls_txt = str(res.get("classification", "")).lower()
    if any(neg in cls_txt for neg in ("benign", "neutral", "polymorphism", "passenger", "sub-threshold", "subthreshold", "non-driver", "compet")):
        return 0
    if any(pos in cls_txt for pos in ("driver", "pathogenic", "loss", "gain", "splice", "polya", "disrupt")):
        return 1
    return 0


def load_engine(server_dir):
    sys.path.insert(0, os.path.abspath(server_dir))
    import server
    def score(wt, mut):
        res = server.run_v33_prediction({"wt_seq": wt, "mut_seq": mut, "flank_bp": FLANK_BP})
        if not res:
            return float("nan"), 0, None
        if "S_raw" in res:
            s = float(res["S_raw"])
        elif "score" in res:
            s = float(res["score"])
        else:
            raise KeyError(f"No S_raw/score in engine output; keys = {sorted(res)}")
        return s, is_frozen_driver_call(res), res
    return score, server


def build_flanks(fa, chrom, pos, ref, alt):   # same as heldout_full_pool_eval.py
    name = fa.contig(chrom)
    if name is None:
        return None, "no_contig"
    start0 = pos - 1 - FLANK_BP
    end0 = pos - 1 + len(ref) + FLANK_BP
    dna = fa.fetch(name, start0, end0)
    if dna is None or len(dna) != end0 - start0:
        return None, "flank_out_of_range"
    if dna[FLANK_BP:FLANK_BP + len(ref)] != ref:
        return None, "ref_mismatch"
    return (dna, dna[:FLANK_BP] + alt + dna[FLANK_BP + len(ref):]), None


# ------------------------------------------------------------ reading files
def open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "rt", encoding="utf-8", errors="replace", newline="")


def norm_chrom(c):
    c = str(c).strip()
    return c[3:] if c.lower().startswith("chr") else c


CAND = {
    "chrom": ["chrom", "chr", "chromosome", "#chrom", "contig"],
    "pos": ["pos", "position", "start", "hg38_pos", "hg38_position"],
    "ref": ["ref", "reference", "ref_allele", "reference_allele"],
    "alt": ["alt", "alternate", "alt_allele", "alternate_allele"],
    "label": ["classification", "splicing_classification", "splice_classification", "validation",
              "label", "class", "category"],
    "variant": ["variant", "hg38", "variant_id", "id", "name"],
}
VARIANT_RE = re.compile(r"(?:chr)?([0-9]{1,2}|X|Y|M|MT)[-:_ ](\d+)[-:_ ]([ACGTN]+)[-:_ >]+([ACGTN]+)", re.I)
LOW_RE = re.compile(r"low.?frequency", re.I)
NEG_RE = re.compile(r"\bnormal\b|non.?splice.?altering", re.I)
POS_RE = re.compile(r"splice.?altering", re.I)


def pick(header, names, override):
    low = {h.lower().lstrip("#"): h for h in header}
    if override:
        if override not in header:
            raise SystemExit(f"Column '{override}' not in header: {header}")
        return override
    for n in names:
        if n.lstrip("#") in low:
            return low[n.lstrip("#")]
    return None


# UCSC Table Browser export of the hg38 SpliceVarDB track (bed-like, no label or ref/alt columns):
# the class is encoded in the colour column, the alleles in the name (transcript-strand HGVS).
UCSC_COLORS = {"219,61,61": ("Splice-altering", 1), "57,135,204": ("Normal", 0),
               "128,82,160": ("Low-frequency", None), "140,140,140": ("Conflicting", None)}
HGVS_SNV = re.compile(r"[cng]\.[-+*0-9_?()]+([ACGT])>([ACGT])\s*$")
COMP = str.maketrans("ACGT", "TGCA")


def load_ucsc(path):
    rows, drops, header = [], Counter(), None
    print(f"UCSC track export detected. Colour map: {UCSC_COLORS}")
    with open_text(path) as f:
        for line in f:
            if not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            if header is None:
                header = [h.lstrip("#") for h in t]
                ic = header.index("color")
                continue
            col = t[ic].strip()
            if col not in UCSC_COLORS:
                drops["unknown_colour:" + col] += 1
                continue
            cname, lab = UCSC_COLORS[col]
            if lab is None:
                drops[cname.lower().replace("-", "_") + "_class"] += 1
                continue
            s, e = int(t[1]), int(t[2])
            m = HGVS_SNV.search(t[3])
            if e - s != 1 or not m:
                drops["not_snv_skipped"] += 1
                continue
            rows.append({"chrom": norm_chrom(t[0]), "pos": e, "hg_ref": m.group(1), "hg_alt": m.group(2), "label": lab})
    return rows, drops


def resolve_alleles(rows, fa, drops):
    """Name alleles are on the transcript strand. Compare the genome base with them: equal -> plus strand,
    equal to the complement -> minus strand (flip both alleles). A base can't match both."""
    out, strand = [], Counter()
    for r in rows:
        name = fa.contig(r["chrom"])
        base = fa.fetch(name, r["pos"] - 1, r["pos"]) if name else None
        hr, ha = r.pop("hg_ref"), r.pop("hg_alt")
        if base is None:
            drops["no_contig_or_out_of_range"] += 1
        elif base == hr:
            r["ref"], r["alt"] = hr, ha
            strand["plus"] += 1
            out.append(r)
        elif base == hr.translate(COMP):
            r["ref"], r["alt"] = hr.translate(COMP), ha.translate(COMP)
            strand["minus"] += 1
            out.append(r)
        else:
            drops["ref_matches_neither_strand"] += 1
    print(f"Strand inferred from the genome: {dict(strand)}; "
          f"ref matched neither strand: {drops['ref_matches_neither_strand']}")
    return out, dict(strand)


def load_svdb(path, a):
    with open_text(path) as f:
        first = f.readline()
    if first.startswith("#chrom") and "chromStart" in first and "color" in first:
        return load_ucsc(path)
    rows, drops = [], Counter()
    with open_text(path) as f:
        lines = [l for l in f if l.strip() and not l.startswith("##")]
    if not lines:
        raise SystemExit("SpliceVarDB file is empty")
    head = lines[0].lstrip("#").rstrip("\r\n")
    delim = "\t" if "\t" in head else ","
    header = next(csv.reader([head], delimiter=delim))
    cols = {k: pick(header, CAND[k], getattr(a, "col_" + k, None)) for k in CAND}
    print(f"SpliceVarDB header: {header}")
    print(f"Columns used: {cols}")
    have_coords = all(cols[k] for k in ("chrom", "pos", "ref", "alt"))
    if not have_coords and not cols["variant"]:
        raise SystemExit("Can't find chrom/pos/ref/alt (or a combined variant column). "
                         "Pass --col-chrom --col-pos --col-ref --col-alt.")
    if not cols["label"]:
        raise SystemExit("Can't find the classification column. Pass --col-label.")
    for rec in csv.DictReader(lines[1:], fieldnames=header, delimiter=delim):
        if have_coords:
            chrom, pos, ref, alt = (rec[cols["chrom"]], rec[cols["pos"]], rec[cols["ref"]], rec[cols["alt"]])
        else:
            m = VARIANT_RE.search(rec.get(cols["variant"], "") or "")
            if not m:
                drops["unparsed_variant"] += 1
                continue
            chrom, pos, ref, alt = m.groups()
        try:
            pos = int(pos)
        except ValueError:
            drops["unparsed_variant"] += 1
            continue
        ref, alt = ref.strip().upper(), alt.strip().upper()
        if "," in alt or not re.fullmatch(r"[ACGTN]+", ref) or not re.fullmatch(r"[ACGTN]+", alt):
            drops["multiallelic_or_symbolic"] += 1
            continue
        txt = rec.get(cols["label"], "") or ""
        if LOW_RE.search(txt):
            drops["low_frequency_class"] += 1
            continue
        if NEG_RE.search(txt):
            lab = 0
        elif POS_RE.search(txt):
            lab = 1
        else:
            drops["unrecognised_label:" + txt[:30]] += 1
            continue
        rows.append({"chrom": norm_chrom(chrom), "pos": pos, "ref": ref, "alt": alt, "label": lab})
    return rows, drops


def read_keys(path):
    """Read (chrom, pos, ref, alt) keys from the calibration file.
    Accepts a standard VCF (chrom, pos, id, ref, alt, ...) AND the variant used here, which has an
    extra leading label column (PATHOGENIC, chrom, pos, id, ref, alt, ...)."""
    keys = set()
    with open_text(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            for off in (0, 1):          # 0 = standard VCF, 1 = extra leading label column
                if len(t) < off + 5:
                    continue
                try:
                    k = (norm_chrom(t[off]), int(t[off + 1]), t[off + 3].strip().upper(), t[off + 4].strip().upper())
                except ValueError:
                    continue
                if re.fullmatch(r"[ACGTN]+", k[2]) and re.fullmatch(r"[ACGTN]+", k[3]):
                    keys.add(k)
                    break
    return keys


# ------------------------------------------------------------ metrics
def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    r = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n)
    return (max(0.0, (c - r) / d), min(1.0, (c + r) / d))


def rank_avg(x):
    o = np.argsort(x, kind="mergesort")
    s = x[o]
    r = np.empty(len(x))
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and s[j + 1] == s[i]:
            j += 1
        r[o[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return r


def metrics(rows):
    n = len(rows)
    if n == 0:
        return {"n": 0}
    y = np.array([r["label"] for r in rows])
    s = np.array([r["score"] for r in rows], float)
    c = np.array([r["call"] for r in rows])
    P, N = int(y.sum()), int((1 - y).sum())
    tp = int(((c == 1) & (y == 1)).sum()); fp = int(((c == 1) & (y == 0)).sum())
    tn = int(((c == 0) & (y == 0)).sum()); fn = int(((c == 0) & (y == 1)).sum())
    div = lambda a, b: a / b if b else float("nan")
    mcc_den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out = {"n": n, "n_pos": P, "n_neg": N, "prevalence": P / n,
           "tp": tp, "fp": fp, "tn": tn, "fn": fn,
           "sensitivity": div(tp, P), "sens_ci95": wilson(tp, P),
           "specificity": div(tn, N), "spec_ci95": wilson(tn, N),
           "ppv": div(tp, tp + fp), "npv": div(tn, tn + fn),
           "accuracy": (tp + tn) / n,
           "mcc": ((tp * tn - fp * fn) / mcc_den) if mcc_den else float("nan")}
    if P and N:
        r = rank_avg(s)
        out["auroc_S_raw"] = float((r[y == 1].sum() - P * (P + 1) / 2) / (P * N))
    return out


def show(name, m):
    if not m.get("n"):
        print(f"\n[{name}] no variants")
        return
    print(f"\n[{name}] n={m['n']} (pos {m['n_pos']}, neg {m['n_neg']}, prevalence {m['prevalence']:.1%})")
    print(f"  sensitivity {m['sensitivity']:.3f}  (95% CI {m['sens_ci95'][0]:.3f}-{m['sens_ci95'][1]:.3f})")
    print(f"  specificity {m['specificity']:.3f}  (95% CI {m['spec_ci95'][0]:.3f}-{m['spec_ci95'][1]:.3f})")
    print(f"  PPV {m['ppv']:.3f}  NPV {m['npv']:.3f}  accuracy {m['accuracy']:.3f}  MCC {m['mcc']:.3f}"
          + (f"  AUROC(S_raw) {m['auroc_S_raw']:.3f}" if "auroc_S_raw" in m else ""))


# ------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--svdb", required=True, help="SpliceVarDB download (tsv/csv/vcf, optionally gzipped)")
    ap.add_argument("--fasta", required=True, help="GRCh38 FASTA file or the folder that holds it")
    ap.add_argument("--calib", default="dev_chr1_11.vcf", help="ClinVar calibration set to exclude")
    ap.add_argument("--server-dir", default=".")
    ap.add_argument("--out", default="eval_splicevardb")
    ap.add_argument("--max-len", type=int, default=15, help="max ref/alt length (ClinVar pool used 15)")
    ap.add_argument("--smoke", type=int, default=0, help="score only N variants, evenly spaced across the file")
    for k in ("chrom", "pos", "ref", "alt", "label"):
        ap.add_argument(f"--col-{k}", default=None)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    score, server = load_engine(a.server_dir)
    expected = getattr(server, "EXPECTED_BUILD_PREFIX", None)
    short = getattr(server, "CORE_SHORT_SHA256", None)
    print(f"Engine build {short}; expected {expected}; verified={short == expected}")
    if expected and short != expected:
        raise SystemExit("Build hash does not match the frozen release. Stopping.")

    rows, drops = load_svdb(a.svdb, a)
    fa = Fasta(resolve_fasta(a.fasta))
    strand_info = None
    if rows and "hg_ref" in rows[0]:
        rows, strand_info = resolve_alleles(rows, fa, drops)
    print(f"\nParsed {len(rows)} SpliceVarDB variants; dropped at parse: {dict(drops)}")

    calib = read_keys(a.calib)
    print(f"Calibration set {a.calib}: {len(calib)} variants (your notes say 955; check this matches)")
    if len(calib) < 10:
        raise SystemExit("Calibration set looks empty or unreadable; overlap removal would be meaningless. Stopping.")

    seen, uniq, dup, conflict = {}, [], 0, 0
    for r in rows:
        k = (r["chrom"], r["pos"], r["ref"], r["alt"])
        if k in seen:
            dup += 1
            conflict += int(seen[k] != r["label"])
            continue
        seen[k] = r["label"]
        uniq.append(r)
    rows = uniq
    long_ = [r for r in rows if len(r["ref"]) > a.max_len or len(r["alt"]) > a.max_len]
    rows = [r for r in rows if len(r["ref"]) <= a.max_len and len(r["alt"]) <= a.max_len]
    for r in rows:
        r["in_calib"] = (r["chrom"], r["pos"], r["ref"], r["alt"]) in calib
    n_overlap = sum(r["in_calib"] for r in rows)
    print(f"Duplicates removed: {dup} (conflicting labels: {conflict}); too long (> {a.max_len}): {len(long_)}")
    print(f"OVERLAP with ClinVar calibration set: {n_overlap} variants (excluded from the headline numbers)")
    if a.smoke:
        step = max(1, len(rows) // a.smoke)
        rows = rows[::step][:a.smoke]

    # scoring, cached
    cache_path = os.path.join(a.out, "scores.tsv")
    cache = {}
    if os.path.exists(cache_path) and not a.smoke:
        with open(cache_path) as f:
            next(f)
            for line in f:
                t = line.rstrip("\n").split("\t")
                cache[(t[0], int(t[1]), t[2], t[3])] = (float(t[4]), int(t[5]), t[6])
        print(f"Loaded {len(cache)} cached scores from {cache_path}")
    status = Counter()
    t0 = time.time()
    for i, r in enumerate(rows):
        k = (r["chrom"], r["pos"], r["ref"], r["alt"])
        if k in cache:
            r["score"], r["call"], r["status"] = cache[k]
        else:
            fl, err = build_flanks(fa, r["chrom"], r["pos"], r["ref"], r["alt"])
            if err:
                r["score"], r["call"], r["status"] = float("nan"), 0, err
            else:
                s, c, res = score(*fl)
                r["score"], r["call"], r["status"] = (s, c, "ok") if res else (float("nan"), 0, "engine_none")
        status[r["status"]] += 1
        if (i + 1) % 500 == 0:
            print(f"  scored {i + 1}/{len(rows)}  ({time.time() - t0:.0f}s)", flush=True)
    print(f"Scoring status: {dict(status)}")
    if not a.smoke:
        with open(cache_path, "w") as f:
            f.write("chrom\tpos\tref\talt\tscore\tcall\tstatus\n")
            for r in rows:
                f.write(f"{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t{r['score']}\t{r['call']}\t{r['status']}\n")

    ok = [r for r in rows if r["status"] == "ok"]
    held = [r for r in ok if not r["in_calib"]]
    c1222 = {str(i) for i in range(12, 23)}
    sets = {
        "HEADLINE: all SpliceVarDB, calibration overlap removed": held,
        "chr12-22 only (the same chromosomes as your locked hold-out), overlap removed":
            [r for r in held if r["chrom"] in c1222],
        "chr1-11 + X, overlap removed (thresholds were tuned on these chromosomes, so not chromosome-held-out)":
            [r for r in held if r["chrom"] not in c1222],
        "reference only: all, overlap NOT removed": ok,
    }
    report = {"svdb_file": os.path.basename(a.svdb), "engine_build": short,
              "parse_drops": dict(drops), "strand_inferred": strand_info, "duplicates": dup, "too_long": len(long_),
              "calibration_size": len(calib), "calibration_overlap": n_overlap,
              "scoring_status": dict(status), "results": {}}
    for name, subset in sets.items():
        m = metrics(subset)
        report["results"][name] = m
        show(name, m)
    with open(os.path.join(a.out, "report.json"), "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nSaved {cache_path} and {os.path.join(a.out, 'report.json')}")


if __name__ == "__main__":
    main()