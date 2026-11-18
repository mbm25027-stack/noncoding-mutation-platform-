# -*- coding: utf-8 -*-
"""
nc_genome_features.py -- gene-annotation helper for the gnomAD experiment (exp6).

Reads a GENCODE GTF once (cached as .npz), then classifies any GRCh38 position into a
NONCODING domain with a distance to the relevant anchor:

    utr       inside an annotated exon but outside every CDS     (5'/3'-UTR exons, incl. polyA region)
    splice    intronic, 1..W bp from the nearest exon edge        (dist = bp to the exon edge; 1-2 = canonical)
    promoter  within [TSS-up, TSS+down] and none of the above     (dist = bp to the nearest TSS)

Anything inside a CDS (any isoform) is dropped, so the analysed variants are noncoding on every
annotated isoform, the same rule as the strictly-noncoding ClinVar pool.
Needs only numpy.
"""
from __future__ import annotations

import bisect
from collections import defaultdict
from pathlib import Path

import numpy as np

from nc_exp_common import log, norm_chrom, open_text, resolve_file


def merge_intervals(iv, gap=0):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + gap + 1:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return out


class IntervalSet:
    """Sorted, merged, 1-based inclusive intervals with O(log n) scalar queries."""

    def __init__(self, merged):
        self.s = [a for a, _ in merged]
        self.e = [b for _, b in merged]

    def contains(self, pos):
        i = bisect.bisect_right(self.s, pos) - 1
        return i >= 0 and pos <= self.e[i]

    def distance(self, pos):
        """0 if inside, else bp to the nearest interval edge (large number if empty)."""
        if not self.s:
            return 10 ** 9
        i = bisect.bisect_right(self.s, pos) - 1
        if i >= 0 and pos <= self.e[i]:
            return 0
        d = 10 ** 9
        if i >= 0:
            d = pos - self.e[i]
        if i + 1 < len(self.s):
            d = min(d, self.s[i + 1] - pos)
        return d


class GenomeFeatures:
    def __init__(self, gtf, chroms=None, cache=None, splice_window=100, prom_up=1000, prom_down=100):
        self.W, self.up, self.down = splice_window, prom_up, prom_down
        raw = self._load(gtf, chroms, cache)
        self.exon, self.cds, self.tss, self.union = {}, {}, {}, {}
        for ch, d in raw.items():
            self.exon[ch] = IntervalSet(merge_intervals(map(tuple, d["exon"])))
            self.cds[ch] = IntervalSet(merge_intervals(map(tuple, d["cds"])))
            self.tss[ch] = sorted(set(d["tss"]))
            dom = [(s - splice_window, e + splice_window) for s, e in merge_intervals(map(tuple, d["exon"]))]
            dom += [(p - prom_up, p + prom_down) if st > 0 else (p - prom_down, p + prom_up)
                    for p, st in zip(d["tss_pos"], d["tss_strand"])]
            self.union[ch] = IntervalSet(merge_intervals(dom))

    # ---- parsing + cache -------------------------------------------------
    def _load(self, gtf, chroms, cache):
        cache = Path(cache) if cache else None
        if cache and cache.exists():
            z = np.load(cache, allow_pickle=False)
            names = sorted({k.split("|")[0] for k in z.files})
            log(f"GTF features from cache {cache}")
            return {ch: dict(exon=z[f"{ch}|exon"].tolist(), cds=z[f"{ch}|cds"].tolist(),
                             tss=z[f"{ch}|tss"].tolist(), tss_pos=z[f"{ch}|tsspos"].tolist(),
                             tss_strand=z[f"{ch}|tssstrand"].tolist()) for ch in names
                    if not chroms or ch in chroms}
        path = resolve_file(gtf, (".gtf", ".gtf.gz"))
        log(f"parsing {path} (one-off; the result is cached) ...")
        raw = defaultdict(lambda: dict(exon=[], cds=[], tss=[], tss_pos=[], tss_strand=[]))
        with open_text(path) as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 8)
                if len(f) < 8 or f[2] not in ("exon", "CDS", "transcript"):
                    continue
                ch = norm_chrom(f[0])
                if chroms and ch not in chroms:
                    continue
                s, e = int(f[3]), int(f[4])
                d = raw[ch]
                if f[2] == "exon":
                    d["exon"].append((s, e))
                elif f[2] == "CDS":
                    d["cds"].append((s, e))
                else:
                    st = 1 if f[6] == "+" else -1
                    p = s if st > 0 else e
                    d["tss"].append(p)
                    d["tss_pos"].append(p)
                    d["tss_strand"].append(st)
        raw = dict(raw)
        if cache:
            arrs = {}
            for ch, d in raw.items():
                arrs[f"{ch}|exon"] = np.array(d["exon"], dtype=np.int64).reshape(-1, 2)
                arrs[f"{ch}|cds"] = np.array(d["cds"], dtype=np.int64).reshape(-1, 2)
                arrs[f"{ch}|tss"] = np.array(d["tss"], dtype=np.int64)
                arrs[f"{ch}|tsspos"] = np.array(d["tss_pos"], dtype=np.int64)
                arrs[f"{ch}|tssstrand"] = np.array(d["tss_strand"], dtype=np.int64)
            np.savez_compressed(cache, **arrs)
        return raw

    # ---- queries -----------------------------------------------------------
    def prefilter(self, chrom, pos):
        u = self.union.get(chrom)
        return bool(u and u.contains(pos))

    def classify(self, chrom, pos):
        """-> (domain, dist) or None (coding or outside every domain)."""
        if chrom not in self.exon or self.cds[chrom].contains(pos):
            return None
        if self.exon[chrom].contains(pos):
            return "utr", 0
        d = self.exon[chrom].distance(pos)
        if 1 <= d <= self.W:
            return "splice", d
        t = self.tss[chrom]
        i = bisect.bisect_left(t, pos)
        cand = [abs(pos - t[j]) for j in (i - 1, i) if 0 <= j < len(t)]
        if cand and min(cand) <= max(self.up, self.down):
            dist = min(cand)
            # honour the strand-specific window used for the prefilter
            return ("promoter", dist) if self.union[chrom].contains(pos) else None
        return None
