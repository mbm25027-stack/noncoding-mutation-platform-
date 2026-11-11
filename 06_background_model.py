"""
06_background_model.py

A context-aware background mutation model in CODING space, with positive and
negative controls. It answers: "which genes and which residues carry MORE
mutations than chance explains?"

The model (a simplified dN/dS-style approach, in the spirit of dNdScv)
  1. Context rates. For each of the 96 trinucleotide classes, the mutation rate
     per site is estimated from ALL coding mutations jointly: synonymous sites
     count fully, missense and nonsense sites count with one genome-wide
     neutral weight each (w_mis, w_non = the overall non-silent : silent ratio,
     estimated from the data; most genes are passengers, so they set it).
     (v2 change. v1 used synonymous mutations only, which overestimates the rate
     of classes that have almost no synonymous sites; see fit_rates.)
  2. Gene rates. Each gene's own background level (replication timing, expression,
     chromatin) is estimated from its synonymous mutations. Gene rates are modelled
     as gamma-distributed (mean 1) with the spread estimated from the data, so a
     gene with few synonymous mutations has an uncertain rate.
  3. Expectation. For any gene or residue, the expected number of mutations =
     number of patients x gene rate x sum of the context rates of the possible
     changes at that gene or residue.
  4. Test. Observed non-silent mutations are compared with the expectation with
     a negative binomial test (Poisson counts + uncertain gene rate), which keeps
     the test honest. False-discovery control is applied over all genes, and over
     ALL codons in the coding genome for the residue test.

Why this step matters for the noncoding project
  If this model rediscovers BRAF V600, NRAS Q61 and other known melanoma hotspots
  WITHOUT flagging huge passenger genes such as TTN and MUC16, the machinery (context
  rates + local rate + proper multiple-testing) works. The same machinery is what a
  noncoding hotspot caller needs on whole-genome data.

Inputs
    results/context/mutations_with_context.csv.gz   from 05_context_and_gene_types.py
    GRCh38.d1.vd1.fa and gencode.v36.annotation.gtf (found automatically)

Outputs
    tables  19_gene_level_test, 20_codon_hotspots, 21_context_rates, 22_control_genes
    figures fig12_qq_gene_level, fig13_observed_vs_expected, fig14_top_hotspots

Limitations (state these in any write-up)
  * Assumes uniform exome coverage; real coverage varies between exons.
  * One transcript per gene (APPRIS principal / CCDS / longest). Residue numbers
    follow that transcript, so they can differ from the MAF's numbering.
  * One dispersion parameter for all genes, and only trinucleotide context + gene-level
    rate. Covariates (replication timing, expression) and finer rate variation
    (dinucleotides, local chromatin) are not modelled.
  * Hypermutated tumours dominate the rate estimates.
Needs: pip install pyfaidx
"""
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # write files, never open windows
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from pyfaidx import Fasta
except ImportError:
    raise SystemExit("pyfaidx is not installed. Run:  python -m pip install pyfaidx  "
                     "(with the same Python you use to run this script).")

# ========================= CHANGE THESE IF NEEDED ===========================
BASE = Path(r"D:\noncoding mutation finding surver")   # <<< project folder
FASTA = None   # <<< None = find the genome FASTA automatically, or give a full path in r"..."
GTF = None     # <<< None = find the GENCODE GTF automatically, or give a full path
ALPHA = None            # <<< gene-rate dispersion: None = estimate from the data; smaller number = more gene-to-gene variation
FDR = 0.05              # <<< significance threshold on the q-value
MIN_PATIENTS_REPORT = 3 # <<< residues need at least this many mutated patients to be listed as hotspots
POSITIVE_CONTROLS = ["BRAF", "NRAS", "TP53", "CDKN2A", "PPP6C", "RAC1", "MAP2K1", "IDH1"]   # <<< expected drivers
NEGATIVE_CONTROLS = ["TTN", "MUC16", "DNAH5", "PCLO", "LRP1B", "CSMD1", "OBSCN", "SYNE1", "RYR2", "FLG"]  # <<< big passenger genes
# ============================================================================

RES = BASE / "results"
TAB = RES / "tables"
FIG = RES / "figures"
LN10 = math.log(10)

# ---- base codes: A=0 C=1 G=2 T=3 (complement of b is 3-b) ----
LUT = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    LUT[ord(_b)] = _i
    LUT[ord(_b.lower())] = _i
BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}

# ---- 96 trinucleotide classes, same order as script 05 ----
SPECTRUM = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CHANNELS_96 = [f"{f}[{s}]{t}" for s in SPECTRUM for f in "ACGT" for t in "ACGT"]
SUBTAB = np.full((4, 4), -1, dtype=np.int16)          # [pyrimidine ref, alt] -> substitution type
SUBTAB[1, 0], SUBTAB[1, 2], SUBTAB[1, 3] = 0, 1, 2    # C>A, C>G, C>T
SUBTAB[3, 0], SUBTAB[3, 1], SUBTAB[3, 2] = 3, 4, 5    # T>A, T>C, T>G

# ---- standard genetic code ----
AA_LETTERS = "ACDEFGHIKLMNPQRSTVWY*"
STOP = AA_LETTERS.index("*")
_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"  # codon order TCAG
AA_IDX = np.zeros(64, dtype=np.int8)
for _a, _b1 in enumerate("TCAG"):
    for _b, _b2 in enumerate("TCAG"):
        for _c, _b3 in enumerate("TCAG"):
            _codon = BASE_CODE[_b1] * 16 + BASE_CODE[_b2] * 4 + BASE_CODE[_b3]
            AA_IDX[_codon] = AA_LETTERS.index(_CODE[_a * 16 + _b * 4 + _c])
SYN, MIS, NON, OTHER = 0, 1, 2, 3
MULT = np.array([16, 4, 1], dtype=np.int16)


# ------------------------------ small helpers -------------------------------
def save_table(df, name, index=True):
    path = TAB / f"{name}.csv"
    df.to_csv(path, index=index)
    print(f"  table   {path.name}")


def save_fig(fig, name):
    path = FIG / f"{name}.png"
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  figure  {path.name}")


def find_file(configured, names, what):
    """Find a real file; Windows extraction often leaves a FOLDER named like the file."""
    if configured:
        p = Path(configured)
        if p.is_file():
            return p
        raise SystemExit(f"{what}: {p} is not a file.")
    for name in names:
        hits = [p for p in BASE.rglob(name) if p.is_file()]
        if hits:
            return sorted(hits, key=lambda p: -p.stat().st_size)[0]
    raise SystemExit(f"Could not find {what} under {BASE}. Looked for: {names}. "
                     "Set FASTA / GTF at the top of the script.")


def nb_nl10_sf(k, shape, rate, lam0):
    """-log10 P(X >= k) when X | g ~ Poisson(g * lam0) and the gene rate g ~ Gamma(shape, rate).
    This is a negative binomial. Computed in log space so tiny p-values do not underflow.

    shape and rate come from the gene's own synonymous mutations (posterior of its rate),
    lam0 is the expected count at a gene rate of exactly 1."""
    k = int(k)
    if k <= 0:
        return 0.0
    lam0 = max(float(lam0), 1e-300)
    p = rate / (rate + lam0)
    q = lam0 / (rate + lam0)
    mean = shape * q / p

    def logpmf(n):
        return (math.lgamma(n + shape) - math.lgamma(shape) - math.lgamma(n + 1)
                + shape * math.log(p) + n * math.log(q))

    if k > mean:
        s = t = 1.0
        for j in range(1, 200000):
            t *= (shape + k + j - 1) / (k + j) * q
            s += t
            if t < 1e-16 * s:
                break
        logp = logpmf(k) + math.log(s)
    else:  # p is large here; use the lower tail
        n = k - 1
        s = t = 1.0
        for j in range(1, n + 1):
            t *= (n - j + 1) / ((n - j + shape) * q)
            s += t
            if t < 1e-16 * s:
                break
        cdf = min(math.exp(logpmf(n)) * s, 1.0)
        logp = math.log(max(1.0 - cdf, 1e-300))
    return float(min(max(-logp / LN10, 0.0), 700.0))


def estimate_alpha(obs, exp):
    """Spread of the gene rate, from how much synonymous counts vary between genes
    beyond Poisson noise: Var(rate) = 1/alpha. Every gene with at least 2 expected
    synonymous mutations counts equally, so one huge gene cannot dominate the estimate.
    Returns alpha clipped to [1, 100]."""
    ok = exp >= 2
    o, e = obs[ok], exp[ok]
    var_g = max(float(np.mean(((o - e) ** 2 - o) / e ** 2)), 0.01)
    return float(np.clip(1.0 / var_g, 1.0, 100.0))


def bh_nlq(nl10_p, m_total=None):
    """Benjamini-Hochberg q-values, returned as -log10(q), computed in log space.
    m_total = total number of tests (untested sites count as p = 1)."""
    nl = np.asarray(nl10_p, dtype=float)
    m = m_total if m_total is not None else len(nl)
    order = np.argsort(-nl, kind="stable")           # most significant first
    ranks = np.arange(1, len(nl) + 1)
    nlq_sorted = nl[order] - math.log10(m) + np.log10(ranks)
    nlq_sorted = np.maximum.accumulate(nlq_sorted[::-1])[::-1]   # q is the running minimum from the tail
    nlq = np.empty_like(nl)
    nlq[order] = np.maximum(nlq_sorted, 0.0)
    return nlq


# ----------------------------- CDS (coding) models --------------------------
def load_cds_models(gtf_path):
    """Read the CDS lines of the GENCODE GTF and pick one transcript per gene.

    Preference: APPRIS principal_1, any APPRIS principal, CCDS tag, then longest
    coding sequence. Only transcripts whose CDS is a whole number of codons and
    starts in frame are kept. Mitochondrial genes are skipped (different genetic code)."""
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    r_tid = re.compile(r'transcript_id "([^"]+)"')
    r_gn = re.compile(r'gene_name "([^"]+)"')
    r_tag = re.compile(r'tag "([^"]+)"')
    tx = {}
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "\tCDS\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS" or f[0] in ("chrM", "chrMT", "MT", "M"):
                continue
            attr = f[8]
            tid, gn = r_tid.search(attr), r_gn.search(attr)
            if not tid or not gn:
                continue
            rec = tx.setdefault(tid.group(1), {"gene": gn.group(1), "chrom": f[0], "strand": f[6],
                                               "segs": [], "tags": set()})
            rec["segs"].append((int(f[3]), int(f[4]), f[7]))
            rec["tags"].update(r_tag.findall(attr))
    best = {}
    for tid, rec in tx.items():
        segs = sorted(rec["segs"])
        length = sum(e - s + 1 for s, e, _ in segs)
        first = segs[0] if rec["strand"] == "+" else segs[-1]
        if length % 3 != 0 or first[2] != "0":
            continue
        tags = rec["tags"]
        score = (int("appris_principal_1" in tags), int(any(t.startswith("appris_principal") for t in tags)),
                 int("CCDS" in tags), length)
        cur = best.get(rec["gene"])
        if cur is None or score > cur["score"]:
            best[rec["gene"]] = {"tid": tid, "chrom": rec["chrom"], "strand": rec["strand"],
                                 "segs": segs, "length": length, "score": score}
    print(f"  coding genes with a usable CDS model: {len(best):,}")
    return best


# ------------------------------ site classification -------------------------
def class_index(ref, alt, left, right):
    """Index (0-95) of the trinucleotide class of a change, folded onto the pyrimidine strand."""
    ref, alt, left, right = (np.asarray(x, dtype=np.int16) for x in (ref, alt, left, right))
    pur = (ref == 0) | (ref == 2)
    r2 = np.where(pur, 3 - ref, ref)
    a2 = np.where(pur, 3 - alt, alt)
    five = np.where(pur, 3 - right, left)
    three = np.where(pur, 3 - left, right)
    return (SUBTAB[r2, a2] * 16 + five * 4 + three).astype(np.int16)


def effect(aa_old, aa_new):
    eff = np.full(aa_old.shape, MIS, dtype=np.int8)
    eff[aa_new == aa_old] = SYN
    eff[(aa_new == STOP) & (aa_old != STOP)] = NON
    eff[(aa_old == STOP) & (aa_new != STOP)] = OTHER
    return eff


# ------------------- enumerate every possible coding change -----------------
def build_opportunities(models, genome, muts):
    """For every coding gene, list all possible single-base changes with their
    trinucleotide class and effect (synonymous / missense / nonsense), and locate
    the observed mutations among them.

    WHY: the model needs to know how many opportunities of each kind a gene offers
    (a big gene offers more), and what each observed mutation does to the protein."""
    by_gene = {g: np.asarray(ix) for g, ix in muts.groupby("Hugo_Symbol").groups.items()}
    POS, REFC, ALTC = muts["pos"].to_numpy(), muts["ref_c"].to_numpy(), muts["alt_c"].to_numpy()
    CHROM, PAT = muts["Chromosome"].to_numpy(), muts["patient"].to_numpy()

    by_chrom = defaultdict(list)
    for g, m in models.items():
        by_chrom[m["chrom"]].append(g)

    opp, info, obs_parts = {}, {}, []
    n_skipped = n_codons = 0
    genome_keys = set(genome.keys())
    for chrom in sorted(by_chrom):
        if chrom not in genome_keys:
            n_skipped += len(by_chrom[chrom])
            continue
        raw = genome[chrom][:].seq.encode("ascii")
        arr = LUT[np.frombuffer(raw, dtype=np.uint8)]
        del raw
        for gene in by_chrom[chrom]:
            m = models[gene]
            gpos = np.concatenate([np.arange(s, e + 1, dtype=np.int64) for s, e, _ in m["segs"]])
            L = len(gpos)
            bg, lg, rg = arr[gpos - 1], arr[gpos - 2], arr[gpos]
            if (bg < 0).any() or (lg < 0).any() or (rg < 0).any():
                n_skipped += 1
                continue
            plus = m["strand"] == "+"
            if plus:
                tb = bg.astype(np.int16)
            else:
                bg, lg, rg = bg[::-1], lg[::-1], rg[::-1]       # transcript order
                tb = (3 - bg).astype(np.int16)
            ci_all = np.arange(L) // 3
            codon_int = (tb[0::3] * 16 + tb[1::3] * 4 + tb[2::3]).astype(np.int16)
            old_codon = codon_int[ci_all]
            mult = MULT[np.arange(L) % 3]
            cls = np.empty((L, 3), dtype=np.int16)
            eff = np.empty((L, 3), dtype=np.int8)
            aa_old = AA_IDX[old_codon]
            for k in (1, 2, 3):
                alt_tx = (tb + k) & 3
                eff[:, k - 1] = effect(aa_old, AA_IDX[old_codon - tb * mult + alt_tx * mult])
                alt_g = alt_tx if plus else 3 - alt_tx
                cls[:, k - 1] = class_index(bg, alt_g, lg, rg)
            opp[gene] = np.bincount((cls.astype(np.int32) * 4 + eff).ravel(), minlength=384) \
                .reshape(96, 4).astype(np.int32)
            n_codons += L // 3
            info[gene] = (m["tid"], L // 3)

            # ---- locate observed mutations of this gene inside the CDS ----
            idx = by_gene.get(gene)
            if idx is None:
                continue
            idx = idx[CHROM[idx] == chrom]
            if len(idx) == 0:
                continue
            gp = gpos if plus else gpos[::-1]
            pos = POS[idx]
            ii = np.minimum(np.searchsorted(gpos, pos), L - 1)
            ok = gpos[ii] == pos
            idx, ii = idx[ok], ii[ok]
            j = ii if plus else L - 1 - ii
            gref = (bg if plus else bg)[j]            # genomic reference base at those sites
            alt_g = ALTC[idx]
            keep = (gref == REFC[idx]) & (alt_g != REFC[idx])
            idx, j, alt_g = idx[keep], j[keep], alt_g[keep]
            if len(idx) == 0:
                continue
            alt_tx = alt_g if plus else 3 - alt_g
            k = (alt_tx - tb[j]) & 3
            ci = j // 3
            new_codon = codon_int[ci] - tb[j] * MULT[j % 3] + alt_tx * MULT[j % 3]
            obs_parts.append({
                "gene": np.full(len(idx), gene, dtype=object), "patient": PAT[idx], "ci": ci,
                "cls": cls[j, k - 1], "eff": eff[j, k - 1],
                "aa_old": AA_IDX[codon_int[ci]], "aa_new": AA_IDX[new_codon],
                "codon_cls": cls.reshape(-1, 3, 3)[ci].reshape(len(idx), 9),
                "codon_eff": eff.reshape(-1, 3, 3)[ci].reshape(len(idx), 9),
            })
        del arr
        print(f"    processed {chrom}", end="\r")
    print(" " * 40, end="\r")
    print(f"  genes modelled: {len(opp):,}   codons: {n_codons:,}   genes skipped: {n_skipped}")

    obs = pd.DataFrame({c: np.concatenate([p[c] for p in obs_parts])
                        for c in ("gene", "patient", "ci", "cls", "eff", "aa_old", "aa_new")})
    CC = np.concatenate([p["codon_cls"] for p in obs_parts])
    CE = np.concatenate([p["codon_eff"] for p in obs_parts])
    return opp, info, obs, CC, CE, n_codons


# ---------------------------------- the model -------------------------------
def fit_rates(opp, obs, n_pat):
    """Per-class mutation rate, estimated jointly from synonymous and non-silent mutations.

    WHY (v2 change): v1 used synonymous mutations only. The genetic code makes some
    trinucleotide classes almost absent at synonymous sites (a class can have 1-3
    synonymous sites but thousands of missense sites). With 0 observed synonymous
    mutations such a class got a huge rate (pseudo-count / tiny number of sites), which
    was then multiplied by thousands of non-silent sites. On simulated data that
    follows the model exactly, expected non-silent mutations came out 1.5x too high.

    Now:   observed_c = rate_c * N * ( syn_sites_c + w_mis * mis_sites_c + w_non * non_sites_c )
    rate_c is the rate at a synonymous site; w_mis and w_non are genome-wide neutral
    dN/dS ratios, estimated by iteration (one number each, not per gene)."""
    S = lambda e: sum(o[:, e] for o in opp.values()).astype(float)
    sites = {SYN: S(SYN), MIS: S(MIS), NON: S(NON)}
    cnt = {e: np.bincount(obs.loc[obs["eff"] == e, "cls"], minlength=96).astype(float) for e in (SYN, MIS, NON)}
    w_mis = w_non = 1.0
    for _ in range(500):
        denom = n_pat * (sites[SYN] + w_mis * sites[MIS] + w_non * sites[NON])
        rate = (cnt[SYN] + cnt[MIS] + cnt[NON] + 0.5) / np.maximum(denom, 1e-9)
        new_mis = cnt[MIS].sum() / (n_pat * (rate * sites[MIS]).sum())
        new_non = cnt[NON].sum() / (n_pat * (rate * sites[NON]).sum())
        done = abs(new_mis - w_mis) + abs(new_non - w_non) < 1e-9
        w_mis, w_non = new_mis, new_non
        if done:
            break
    rate[denom <= 0] = np.median(rate[denom > 0])
    t = pd.DataFrame({"synonymous_sites": sites[SYN].astype(int), "synonymous_mutations": cnt[SYN].astype(int),
                      "missense_sites": sites[MIS].astype(int), "missense_mutations": cnt[MIS].astype(int),
                      "nonsense_sites": sites[NON].astype(int), "nonsense_mutations": cnt[NON].astype(int),
                      "rate_per_site_per_patient": rate}, index=CHANNELS_96)
    print(f"  synonymous mutations: {int(cnt[SYN].sum()):,} "
          f"(fewest in one class: {int(cnt[SYN].min())}, median {int(np.median(cnt[SYN]))})")
    print(f"  neutral non-silent : silent weights (genome-wide)   missense w = {w_mis:.3f}   nonsense w = {w_non:.3f}")
    print("     (about 1 means the trinucleotide model matches the data; clearly below 1 means fewer "
          "non-silent mutations than silent ones predict, e.g. purifying selection)")
    return rate, t, w_mis, w_non


def gene_level(opp, info, obs, rate, n_pat, w_mis, w_non):
    """Gene table: observed versus expected non-silent mutations.

    WHY: a mutation count per gene is useless on its own (huge genes win). This
    table adds the expectation from sequence context and the gene's own background
    level (with its uncertainty), so only genes with a real excess stand out."""
    genes = list(opp)
    O = np.stack([opp[g] for g in genes]).astype(float)           # genes x 96 x 4
    exp_syn = n_pat * (O[:, :, SYN] @ rate)
    exp_ns_base = n_pat * ((w_mis * O[:, :, MIS] + w_non * O[:, :, NON]) @ rate)
    obs_syn = obs[obs["eff"] == SYN].groupby("gene").size().reindex(genes).fillna(0).to_numpy()
    obs_ns = obs[obs["eff"].isin([MIS, NON])].groupby("gene").size().reindex(genes).fillna(0).to_numpy()
    alpha = ALPHA if ALPHA is not None else estimate_alpha(obs_syn, exp_syn)
    print(f"  gene-rate dispersion alpha = {alpha:.1f} "
          f"(gene rates vary by about +/-{100 / math.sqrt(alpha):.0f}% beyond chance)")
    post_shape = alpha + obs_syn
    post_rate = alpha + exp_syn
    g_mult = post_shape / post_rate
    nl10 = np.array([nb_nl10_sf(o, a, b, l) for o, a, b, l in zip(obs_ns, post_shape, post_rate, exp_ns_base)])
    exp_ns = g_mult * exp_ns_base
    t = pd.DataFrame({"transcript": [info[g][0] for g in genes], "codons": [info[g][1] for g in genes],
                      "gene_rate_multiplier": g_mult.round(3),
                      "obs_synonymous": obs_syn.astype(int), "exp_synonymous": exp_syn.round(2),
                      "obs_nonsilent": obs_ns.astype(int), "exp_nonsilent": exp_ns.round(2),
                      "obs_over_exp": (obs_ns / np.maximum(exp_ns, 1e-9)).round(2),
                      "nl10_p": nl10.round(3), "nl10_q": bh_nlq(nl10).round(3),
                      "post_shape": post_shape, "post_rate": post_rate}, index=pd.Index(genes, name="gene"))
    t["q_value"] = 10.0 ** (-t["nl10_q"])
    return t.sort_values("nl10_p", ascending=False)


def codon_level(obs, CC, CE, rate, n_pat, gt, n_codons_total, w_mis, w_non):
    """Residue (codon) table: observed versus expected non-silent mutations.

    WHY: this is hotspot calling in coding space. Expected = patients x gene rate x
    sum of context rates of the possible amino-acid-changing substitutions of that
    codon. Multiple testing is applied over EVERY codon of the coding genome, so a
    residue needs a very strong excess to pass."""
    pos = np.flatnonzero(obs["eff"].isin([MIS, NON]).to_numpy())
    ns = obs.iloc[pos].copy()
    ns["pos"] = pos
    ns["change"] = [f"{AA_LETTERS[a]}{c + 1}{AA_LETTERS[b]}"
                    for a, c, b in zip(ns["aa_old"], ns["ci"], ns["aa_new"])]
    grp = ns.groupby(["gene", "ci"])
    t = grp.agg(events=("pos", "size"), first=("pos", "first"), aa_old=("aa_old", "first")).reset_index()
    w_eff = np.array([0.0, w_mis, w_non, 0.0])                # weight by effect: SYN, MIS, NON, OTHER
    first = t["first"].to_numpy()
    unit = (w_eff[CE[first]] * rate[CC[first]]).sum(axis=1)
    lam0 = n_pat * unit                                       # expected at a gene rate of exactly 1
    shape, rate_g = t["gene"].map(gt["post_shape"]).to_numpy(), t["gene"].map(gt["post_rate"]).to_numpy()
    t["expected"] = lam0 * shape / rate_g
    t["nl10_p"] = [nb_nl10_sf(k, a, b, l) for k, a, b, l in zip(t["events"], shape, rate_g, lam0)]
    t["nl10_q"] = bh_nlq(t["nl10_p"].to_numpy(), m_total=n_codons_total)
    t["residue"] = [f"{AA_LETTERS[a]}{c + 1}" for a, c in zip(t["aa_old"], t["ci"])]
    t["obs_over_exp"] = (t["events"] / np.maximum(t["expected"], 1e-9)).round(1)
    seen = (ns.groupby(["gene", "ci", "change"]).size().reset_index(name="n")
              .sort_values(["gene", "ci", "n"], ascending=[True, True, False]))
    seen["txt"] = seen["change"] + " x" + seen["n"].astype(str)
    names = seen.groupby(["gene", "ci"])["txt"].apply(lambda s: "; ".join(s.head(4)))
    t = t.join(names.rename("changes_seen"), on=["gene", "ci"])
    t["q_value"] = 10.0 ** (-t["nl10_q"])
    t["expected"] = t["expected"].round(4)
    t = t.rename(columns={"ci": "codon_index0"}).drop(columns=["first", "aa_old"])
    return t.sort_values(["nl10_q", "events"], ascending=False)[
        ["gene", "residue", "events", "expected", "obs_over_exp", "nl10_p", "nl10_q", "q_value", "changes_seen",
         "codon_index0"]]


# ---------------------------------- figures ---------------------------------
def fig_qq(gt):
    """Figure 12. QQ plot of gene-level p-values.
    WHY: a well-calibrated null model hugs the diagonal; the few points that lift
    off are the genes with a real excess. A curve lifting off everywhere would mean
    the model is miscalibrated. Poisson p-values are discrete, so some conservatism is normal."""
    p = np.sort(gt["nl10_p"].to_numpy())[::-1]
    n = len(p)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(expected, p, s=6, color="#1F3864")
    lim = max(expected.max(), p.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="#E76F51", lw=1)
    for g in POSITIVE_CONTROLS:
        if g in gt.index:
            ax.annotate(g, (expected[np.searchsorted(-p, -gt.loc[g, "nl10_p"])], gt.loc[g, "nl10_p"]),
                        fontsize=7, color="#E76F51", xytext=(4, 0), textcoords="offset points")
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title("Gene-level test: observed vs expected p-values", fontweight="bold", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_fig(fig, "fig12_qq_gene_level")


def fig_obs_exp(gt):
    """Figure 13. Observed versus expected non-silent mutations per gene.
    WHY: the central picture. Big passenger genes (TTN, MUC16) have many mutations,
    but they sit on the diagonal once expectation is modelled. Real drivers sit far above it."""
    d = gt[(gt["exp_nonsilent"] > 0)]
    fig, ax = plt.subplots(figsize=(6.5, 6))
    ax.scatter(d["exp_nonsilent"].clip(lower=1e-2), d["obs_nonsilent"].clip(lower=0.5), s=5, color="#8D99AE", alpha=0.6)
    lim = [max(min(d["exp_nonsilent"].min(), 1.0) * 0.5, 1e-2),
           max(d["exp_nonsilent"].max(), d["obs_nonsilent"].max()) * 2]
    ax.plot(lim, lim, color="#444444", lw=1, ls="--")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    for names, color in [(POSITIVE_CONTROLS, "#E76F51"), (NEGATIVE_CONTROLS, "#1F3864")]:
        for g in names:
            if g in gt.index:
                r = gt.loc[g]
                ax.scatter(max(r["exp_nonsilent"], 1e-2), max(r["obs_nonsilent"], 0.5), s=22, color=color, zorder=3)
                ax.annotate(g, (max(r["exp_nonsilent"], 1e-2), max(r["obs_nonsilent"], 0.5)), fontsize=7,
                            xytext=(3, 3), textcoords="offset points", color=color)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Expected non-silent mutations (model)")
    ax.set_ylabel("Observed non-silent mutations")
    ax.set_title("Observed vs expected per gene\n(orange = expected drivers, blue = large passenger genes)",
                 fontweight="bold", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_fig(fig, "fig13_observed_vs_expected")


def fig_hotspots(ct, top=15):
    """Figure 14. Top residues by significance.
    WHY: shows which residues the model calls hotspots and how far above
    expectation they lie."""
    d = ct.head(top).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.barh([f"{g} {r}" for g, r in zip(d["gene"], d["residue"])], d["events"], color="#E76F51")
    for y, (n, e) in enumerate(zip(d["events"], d["expected"])):
        ax.text(n + 0.5, y, f"{n} obs vs {e:.3f} exp", va="center", fontsize=8)
    ax.set_xlim(0, d["events"].max() * 1.35)
    ax.set_xlabel("Mutations at the residue (patients)")
    ax.set_title("Most significant hotspot residues", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save_fig(fig, "fig14_top_hotspots")


# ----------------------------------- main -----------------------------------
def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    src = RES / "context" / "mutations_with_context.csv.gz"
    if not src.exists():
        raise SystemExit(f"Missing {src}. Run 05_context_and_gene_types.py first.")
    fasta_path = find_file(FASTA, ["GRCh38.d1.vd1.fa", "GRCh38*.fa", "*.fa", "*.fasta"], "reference FASTA")
    gtf_path = find_file(GTF, ["gencode.v36.annotation.gtf", "gencode*.gtf", "gencode*.gtf.gz"], "GENCODE GTF")
    print(f"FASTA: {fasta_path}\nGTF:   {gtf_path}\n")

    cols = ["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
            "Variant_Type", "patient"]
    m = pd.read_csv(src, usecols=cols, low_memory=False)
    n_pat = m["patient"].nunique()
    m = m[(m["Variant_Type"] == "SNP")].copy()
    m["ref_c"] = m["Reference_Allele"].map(BASE_CODE)
    m["alt_c"] = m["Tumor_Seq_Allele2"].map(BASE_CODE)
    m = m.dropna(subset=["ref_c", "alt_c", "Hugo_Symbol"]).reset_index(drop=True)
    m["pos"] = m["Start_Position"].astype(np.int64)
    m["ref_c"] = m["ref_c"].astype(np.int16)
    m["alt_c"] = m["alt_c"].astype(np.int16)
    m["Chromosome"] = m["Chromosome"].astype(str)
    print(f"Patients: {n_pat}   point mutations: {len(m):,}")

    print("\nReading the GENCODE coding sequences (about a minute)")
    models = load_cds_models(gtf_path)

    print("\nListing every possible coding change and placing the observed mutations (a few minutes)")
    genome = Fasta(str(fasta_path))
    opp, info, obs, CC, CE, n_codons = build_opportunities(models, genome, m)
    print(f"  observed point mutations inside the chosen CDS: {len(obs):,}")

    print("\nFitting the background model")
    rate, rate_t, w_mis, w_non = fit_rates(opp, obs, n_pat)
    save_table(rate_t, "21_context_rates")

    print("\nTests")
    gt = gene_level(opp, info, obs, rate, n_pat, w_mis, w_non)
    save_table(gt.drop(columns=["post_shape", "post_rate"]), "19_gene_level_test")
    ct = codon_level(obs, CC, CE, rate, n_pat, gt, n_codons, w_mis, w_non)
    ct_report = ct[ct["events"] >= MIN_PATIENTS_REPORT]
    save_table(ct_report, "20_codon_hotspots", index=False)

    tot_o, tot_e = gt["obs_nonsilent"].sum(), gt["exp_nonsilent"].sum()
    print(f"\nCalibration check: non-silent mutations observed {int(tot_o):,}, expected by the model {tot_e:,.0f} "
          f"(ratio {tot_o / tot_e:.2f}; should be close to 1)")
    sig_g = gt[gt["q_value"] < FDR]
    sig_c = ct_report[ct_report["q_value"] < FDR]
    print(f"\nGenes with q < {FDR}: {len(sig_g)} of {len(gt):,}")
    print(f"Nominal p < 0.05 in {100 * (gt['nl10_p'] > -math.log10(0.05)).mean():.1f}% of genes "
          "(about 5% or less is expected for a well-calibrated test; real drivers add a few)")
    print(f"Residues with q < {FDR} and >= {MIN_PATIENTS_REPORT} patients: {len(sig_c)}")

    print("\nPositive controls (expected drivers): gene-level result and best residue")
    rows = []
    for g in POSITIVE_CONTROLS:
        if g not in gt.index:
            print(f"  {g}: not in the model")
            continue
        r = gt.loc[g]
        c = ct[ct["gene"] == g].head(1)
        best = f"{c['residue'].iloc[0]} ({int(c['events'].iloc[0])} pts, q={c['q_value'].iloc[0]:.1e})" if len(c) else "none"
        print(f"  {g:8s} obs {int(r['obs_nonsilent']):4d}  exp {r['exp_nonsilent']:8.1f}  "
              f"gene q={r['q_value']:.1e}  best residue: {best}")
        rows.append({"role": "positive", "gene": g, "obs": int(r["obs_nonsilent"]), "exp": r["exp_nonsilent"],
                     "gene_q": r["q_value"], "best_residue": best})
    print("\nNegative controls (large passenger genes): should NOT be flagged")
    for g in NEGATIVE_CONTROLS:
        if g not in gt.index:
            continue
        r = gt.loc[g]
        flag = "FLAGGED" if r["q_value"] < FDR else "not flagged"
        print(f"  {g:8s} obs {int(r['obs_nonsilent']):4d}  exp {r['exp_nonsilent']:8.1f}  "
              f"rate x{r['gene_rate_multiplier']:.2f}  gene q={r['q_value']:.2f}  {flag}")
        rows.append({"role": "negative", "gene": g, "obs": int(r["obs_nonsilent"]), "exp": r["exp_nonsilent"],
                     "gene_q": r["q_value"], "best_residue": ""})
    save_table(pd.DataFrame(rows), "22_control_genes", index=False)

    print("\nTop 15 hotspot residues")
    print(ct_report.head(15)[["gene", "residue", "events", "expected", "q_value", "changes_seen"]].to_string(index=False))

    print("\nFigures")
    fig_qq(gt)
    fig_obs_exp(gt)
    fig_hotspots(ct_report if len(ct_report) else ct)
    print("\nDone.")


if __name__ == "__main__":
    main()