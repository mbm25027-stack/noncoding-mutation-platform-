#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nc_exp_common.py -- shared code for the NC-MutScan v3.3 follow-up experiments
(exp1_full_pool_eval.py and exp2_distance_control.py).  Needs only numpy.

Utility command (finds the function to call in your own code):
    python nc_exp_common.py probe server        # lists callables in server.py
"""
from __future__ import annotations

import csv
import gzip
import importlib
import inspect
import math
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# generic helpers
# ---------------------------------------------------------------------------
_trapz = getattr(np, "trapezoid", None) or np.trapz


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def open_text(path):
    """Open plain or gzip/bgzip text (.gz, .bgz or anything else) by checking the file's magic bytes."""
    path = str(path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "rt", encoding="utf-8", errors="replace", newline="")


def resolve_file(path, suffixes):
    """Return a file path; if `path` is a folder (as 'GRCh38.d1.vd1.fa' and
    'gencode.v36.annotation.gtf' are in your project) pick the largest file
    inside it that ends with one of `suffixes`."""
    p = Path(path)
    if p.is_file():
        return p
    if p.is_dir():
        cands = [f for f in p.rglob("*") if f.is_file() and str(f).lower().endswith(tuple(suffixes))]
        if cands:
            best = max(cands, key=lambda f: f.stat().st_size)
            log(f"{path} is a folder; using {best}")
            return best
    raise FileNotFoundError(f"Cannot find a {suffixes} file at {path}")


def norm_chrom(c) -> str:
    c = str(c)
    return c if c.startswith("chr") else "chr" + c


def chrom_set(spec: str | None):
    """'12-22' | '1-11,X' | 'all' -> set of 'chrN' names (None = all)."""
    if not spec or spec.lower() == "all":
        return None
    out = set()
    for tok in spec.split(","):
        tok = tok.strip()
        if "-" in tok and tok.replace("-", "").replace("chr", "").isdigit():
            a, b = tok.replace("chr", "").split("-")
            out.update(f"chr{i}" for i in range(int(a), int(b) + 1))
        else:
            out.add(norm_chrom(tok))
    return out


def chrom_order(c: str) -> int:
    s = c.replace("chr", "")
    return int(s) if s.isdigit() else {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(s, 99)


def write_csv(path, rows, fieldnames=None):
    rows = list(rows)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------
# variant pools (ClinVar VCF or a plain table)
# ---------------------------------------------------------------------------
PATH_TOK = {"pathogenic", "likely_pathogenic"}
BEN_TOK = {"benign", "likely_benign"}
# consequence terms allowed in a "strictly noncoding" record (any other term, e.g.
# missense/synonymous/frameshift/stop/start, removes the record -> isoform purge)
NONCODING_OK = {
    "intron_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "splice_donor_variant",
    "splice_acceptor_variant", "splice_region_variant", "upstream_gene_variant",
    "downstream_gene_variant", "non_coding_transcript_variant",
    "non_coding_transcript_exon_variant", "genic_upstream_transcript_variant",
    "genic_downstream_transcript_variant", "intergenic_variant", "no_sequence_alteration",
    "splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
}
_RX = {k: re.compile(rf"(?:^|;){k}=([^;]+)") for k in ("CLNSIG", "MC", "GENEINFO", "CLNREVSTAT")}


def label_from_clnsig(s: str):
    low = s.lower()
    if "conflicting" in low or "uncertain" in low:
        return None
    toks = {t for t in re.split(r"[|/,]", low) if t}
    p, b = bool(toks & PATH_TOK), bool(toks & BEN_TOK)
    if p and not b:
        return 1
    if b and not p:
        return 0
    return None


def _info(info: str, key: str) -> str:
    m = _RX[key].search(info)
    return m.group(1) if m else ""


def parse_vcf_line(line: str, require_label=True):
    f = line.rstrip("\n").split("\t")
    if len(f) < 8:
        return None
    chrom, pos, vid, ref, alt, info = f[0], f[1], f[2], f[3], f[4], f[7]
    if "," in alt or alt in (".", "*") or ref == ".":
        return None
    lab = label_from_clnsig(_info(info, "CLNSIG"))
    if lab is None and require_label:
        return None
    gene = _info(info, "GENEINFO").split("|")[0].split(":")[0]
    return dict(chrom=norm_chrom(chrom), pos=int(pos), ref=ref.upper(), alt=alt.upper(),
                label=lab, gene=gene, mc=_info(info, "MC"), revstat=_info(info, "CLNREVSTAT"))


def mc_is_noncoding(mc: str) -> bool:
    if not mc:
        return False
    terms = [t.split("|")[-1] for t in mc.split(",") if t]
    return bool(terms) and all(t in NONCODING_OK for t in terms)


def make_id(r) -> str:
    return f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}"


def load_pool(path, chroms=None, max_len=None):
    """Load variants from a ClinVar-style VCF(.gz) (label from CLNSIG) or from a
    TSV/CSV with columns chrom,pos,ref,alt,label[,gene]."""
    path = str(path)
    rows = []
    low = path.lower()
    if low.endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
        with open_text(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                r = parse_vcf_line(line)
                if r:
                    rows.append(r)
    else:
        with open_text(path) as fh:
            sample = fh.read(4096)
            fh.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,;") if sample else csv.excel_tab
            for d in csv.DictReader(fh, dialect=dialect):
                d = {k.strip().lower(): v for k, v in d.items()}
                rows.append(dict(chrom=norm_chrom(d.get("chrom") or d.get("chr")), pos=int(d["pos"]),
                                 ref=d["ref"].upper(), alt=d["alt"].upper(), label=int(float(d["label"])),
                                 gene=d.get("gene", ""), mc=d.get("mc", ""), revstat=d.get("revstat", "")))
    keep = []
    for r in rows:
        if chroms and r["chrom"] not in chroms:
            continue
        if max_len and max(len(r["ref"]), len(r["alt"])) > max_len:
            continue
        r["id"] = make_id(r)
        keep.append(r)
    return keep


def build_pool_from_clinvar(clinvar, chroms, max_len=15, min_stars=0):
    """APPROXIMATE re-creation of the strictly-noncoding ClinVar pool (see README text
    in exp1).  Prefer the pool file your own build_holdout.py produced."""
    stars_map = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
                 "criteria_provided,_multiple_submitters,_no_conflicts": 2,
                 "criteria_provided,_single_submitter": 1, "criteria_provided,_conflicting_classifications": 1,
                 "criteria_provided,_conflicting_interpretations": 1}
    out, seen = [], set()
    n = 0
    with open_text(clinvar) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n += 1
            if n % 1_000_000 == 0:
                log(f"  scanned {n:,} ClinVar records, kept {len(out):,}")
            c = line.split("\t", 1)[0]
            if chroms and norm_chrom(c) not in chroms:
                continue
            r = parse_vcf_line(line)
            if not r or max(len(r["ref"]), len(r["alt"])) > max_len or not mc_is_noncoding(r["mc"]):
                continue
            if min_stars and stars_map.get(r["revstat"].lower(), 0) < min_stars:
                continue
            r["id"] = make_id(r)
            if r["id"] in seen:
                continue
            seen.add(r["id"])
            out.append(r)
    return out


# ---------------------------------------------------------------------------
# FASTA (indexed random access, no dependencies)
# ---------------------------------------------------------------------------
class Fasta:
    def __init__(self, path):
        self.path = resolve_file(path, (".fa", ".fasta", ".fna"))
        self.fh = open(self.path, "rb")
        fai = Path(str(self.path) + ".fai")
        if not fai.exists():
            log(f"building FASTA index {fai} (one-off, a few minutes for a whole genome) ...")
            self._build_fai(fai)
        self.idx = {}
        for line in open(fai, encoding="utf-8"):
            n, ln, off, lb, lw = line.rstrip("\n").split("\t")[:5]
            self.idx[n] = (int(ln), int(off), int(lb), int(lw))

    def _build_fai(self, fai):
        out, name, length, off, lb, lw, pos = [], None, 0, 0, 0, 0, 0
        with open(self.path, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    if name is not None:
                        out.append((name, length, off, lb, lw))
                    name = line[1:].split()[0].decode()
                    length, off, lb, lw = 0, pos + len(line), 0, 0
                else:
                    if lb == 0:
                        lb, lw = len(line.rstrip(b"\r\n")), len(line)
                    length += len(line.rstrip(b"\r\n"))
                pos += len(line)
        if name is not None:
            out.append((name, length, off, lb, lw))
        with open(fai, "w", encoding="utf-8") as fh:
            for r in out:
                fh.write("\t".join(map(str, r)) + "\n")

    def has(self, chrom):
        return chrom in self.idx or chrom.replace("chr", "") in self.idx

    def fetch(self, chrom, start, end):
        """0-based half-open [start, end)."""
        key = chrom if chrom in self.idx else chrom.replace("chr", "")
        ln, off, lb, lw = self.idx[key]
        start, end = max(0, start), min(end, ln)
        if end <= start:
            return ""
        b0 = off + (start // lb) * lw + start % lb
        b1 = off + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        self.fh.seek(b0)
        return self.fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def make_window(fa: Fasta, chrom, pos, ref, alt, flank=25):
    """WT / MUT sequences with `flank` bases either side of the variant.
    Returns (wt, mut, center_index) or None if the reference allele does not match."""
    s0 = pos - 1
    left = fa.fetch(chrom, s0 - flank, s0)
    if fa.fetch(chrom, s0, s0 + len(ref)) != ref:
        return None
    right = fa.fetch(chrom, s0 + len(ref), s0 + len(ref) + flank)
    return left + ref + right, left + alt + right, len(left)


# ---------------------------------------------------------------------------
# engine adapter
# ---------------------------------------------------------------------------
def load_engine(spec: str):
    """spec = 'module:function', default 'my_engine:score' (copy my_engine_template.py)."""
    mod, _, fn = spec.partition(":")
    sys.path.insert(0, os.getcwd())
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        f = getattr(importlib.import_module(mod), fn or "score")
    except ModuleNotFoundError as e:
        raise SystemExit(f"Engine module '{mod}' not found ({e}). Copy my_engine_template.py to "
                         f"{mod}.py and wire it to NC-MutScan v3.3 (see the comments inside).")
    return f


def normalise_result(res) -> dict:
    if not isinstance(res, dict):
        raise TypeError("engine score() must return a dict with s_raw, s_unified, call")
    s_raw = float(res.get("s_raw", res.get("S_raw", 0.0)) or 0.0)
    s_uni = float(res.get("s_unified", res.get("S_unified", 0.0)) or 0.0)
    call = res.get("call")
    call = (s_uni >= 4.2) if call is None else bool(call)
    feats = {f"feat_{k}": float(v) for k, v in (res.get("features") or {}).items()}
    return dict(s_raw=s_raw, s_unified=s_uni, call=int(call), module=str(res.get("module", "")), **feats)


def probe(module_name: str):
    sys.path.insert(0, os.getcwd())
    m = importlib.import_module(module_name)
    print(f"Callables in {module_name}:")
    for n, f in inspect.getmembers(m, inspect.isfunction):
        if f.__module__ == m.__name__:
            doc = (inspect.getdoc(f) or "").splitlines()
            print(f"  {n}{inspect.signature(f)}   {doc[0] if doc else ''}")


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def lr_plus(tp, fn, fp, tn, z=1.96):
    """LR+ with Katz log-method CI (0.5 added to every cell if any cell is zero)."""
    if 0 in (tp, fn, fp, tn):
        tp, fn, fp, tn = tp + .5, fn + .5, fp + .5, tn + .5
    sens, fpr = tp / (tp + fn), fp / (fp + tn)
    lr = sens / fpr
    se = math.sqrt(1 / tp - 1 / (tp + fn) + 1 / fp - 1 / (fp + tn))
    return lr, math.exp(math.log(lr) - z * se), math.exp(math.log(lr) + z * se)


def ppv_at_prior(sens, spec, prior):
    return sens * prior / (sens * prior + (1 - spec) * (1 - prior))


def confusion(y, call, w=None):
    y, call = np.asarray(y, float), np.asarray(call, float)
    w = np.ones_like(y) if w is None else w
    return (float((w * y * call).sum()), float((w * (1 - y) * call).sum()),
            float((w * (1 - y) * (1 - call)).sum()), float((w * y * (1 - call)).sum()))  # tp fp tn fn


class RankPrep:
    """Pre-sorted scores so weighted ROC-AUC / average precision are O(n) per call
    (used for the gene-cluster bootstrap).  Ties are handled exactly."""

    def __init__(self, score):
        score = np.nan_to_num(np.asarray(score, float), nan=-1e18)
        self.order = np.argsort(-score, kind="mergesort")
        s = score[self.order]
        self.bounds = np.r_[np.where(np.diff(s) != 0)[0], len(s) - 1]

    def metrics(self, y, w=None):
        ys = np.asarray(y, float)[self.order]
        ws = np.ones_like(ys) if w is None else np.asarray(w, float)[self.order]
        tp = np.cumsum(ws * ys)[self.bounds]
        fp = np.cumsum(ws * (1 - ys))[self.bounds]
        P, N = tp[-1], fp[-1]
        if P <= 0 or N <= 0:
            return float("nan"), float("nan")
        auc = float(_trapz(np.r_[0, tp / P], np.r_[0, fp / N]))
        prec = tp / np.maximum(tp + fp, 1e-12)
        ap = float(np.sum(np.diff(np.r_[0, tp / P]) * prec))
        return auc, ap


def sens_at_spec(score, y, target_spec):
    score, y = np.asarray(score, float), np.asarray(y)
    neg = np.sort(score[y == 0])
    thr = np.quantile(neg, target_spec) if len(neg) else np.inf
    return float(np.mean(score[y == 1] > thr)) if (y == 1).any() else float("nan")


def cluster_index(keys):
    uniq, inv = np.unique(np.asarray(keys), return_inverse=True)
    return inv, len(uniq)


def cluster_bootstrap(cl_idx, n_cl, B, rng, fn):
    """fn(weights) -> tuple of statistics; weights = how often each variant's gene
    cluster was drawn.  Returns array (B, n_stats)."""
    out = []
    for _ in range(B):
        counts = np.bincount(rng.integers(0, n_cl, n_cl), minlength=n_cl)
        out.append(fn(counts[cl_idx].astype(float)))
    return np.array(out, float)


def pct_ci(arr, q=(2.5, 97.5)):
    return np.nanpercentile(arr, q, axis=0)


def fmt_pct(x, d=2):
    return f"{100 * x:.{d}f}%"


def fmt_ci(est, lo, hi, d=2):
    return f"{100 * est:.{d}f}% ({100 * lo:.{d}f}\u2013{100 * hi:.{d}f})"


def chi2_sf(x, df):
    try:
        from scipy.stats import chi2
        return float(chi2.sf(x, df))
    except Exception:
        a, xx = df / 2.0, x / 2.0            # regularised upper incomplete gamma Q(a, x)
        if xx <= 0:
            return 1.0
        if xx < a + 1:
            s, t, n = 1.0 / a, 1.0 / a, 1
            while abs(t) > abs(s) * 1e-14 and n < 1000:
                t *= xx / (a + n)
                s += t
                n += 1
            return max(0.0, 1.0 - s * math.exp(-xx + a * math.log(xx) - math.lgamma(a)))
        b, c, d = xx + 1 - a, 1e300, 1.0 / (xx + 1 - a)
        h = d
        for i in range(1, 1000):
            an = -i * (i - a)
            b += 2
            d = 1.0 / (an * d + b) if abs(an * d + b) > 1e-300 else 1e300
            c = b + an / c
            h *= d * c
            if abs(d * c - 1) < 1e-14:
                break
        return float(math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * h)


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "probe":
        probe(sys.argv[2])
    else:
        print(__doc__)
