"""
09_annotate_and_score_hotspots.py (Fixed Row Index Alignment)

1. Calibrates window expectation using local regional mutation rates (100kb genomic bins).
2. Annotates noncoding hotspot windows with overlapping GENCODE genes using exact row indices.
3. Produces a prioritized noncoding candidate hotspot table.

Outputs:
  - results/tables/25_prioritized_noncoding_hotspots.csv
"""

import math
from pathlib import Path
import numpy as np
import pandas as pd

# ========================= CONFIGURATION ===========================
BASE = Path(r"D:\noncoding mutation finding surver")
RES = BASE / "results"
TAB = RES / "tables"

LN10 = math.log(10)
# ===================================================================


def find_gtf():
    hits = list(BASE.rglob("gencode*.gtf*"))
    if hits:
        return sorted(hits, key=lambda p: -p.stat().st_size)[0]
    return None


def nb_nl10_sf(k, shape, rate, lam0):
    k = int(k)
    if k <= 0:
        return 0.0
    lam0 = max(float(lam0), 1e-300)
    p = rate / (rate + lam0)
    q = lam0 / (rate + lam0)
    mean = shape * q / p

    def logpmf(n):
        return (
            math.lgamma(n + shape)
            - math.lgamma(shape)
            - math.lgamma(n + 1)
            + shape * math.log(p)
            + n * math.log(q)
        )

    if k > mean:
        s = t = 1.0
        for j in range(1, 200000):
            t *= (shape + k + j - 1) / (k + j) * q
            s += t
            if t < 1e-16 * s:
                break
        logp = logpmf(k) + math.log(s)
    else:
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


def bh_nlq(nl10_p, m_total=None):
    nl = np.asarray(nl10_p, dtype=float)
    m = m_total if m_total is not None else len(nl)
    if m == 0 or len(nl) == 0:
        return np.array([])
    order = np.argsort(-nl, kind="stable")
    ranks = np.arange(1, len(nl) + 1)
    nlq_sorted = nl[order] - math.log10(m) + np.log10(ranks)
    nlq_sorted = np.maximum.accumulate(nlq_sorted[::-1])[::-1]
    nlq = np.empty_like(nl)
    nlq[order] = np.maximum(nlq_sorted, 0.0)
    return nlq


# 1. Load Data
src_win = TAB / "23_noncoding_window_hotspots.csv"
src_snv = RES / "pcawg_melanoma_wgs_snvs.csv.gz"

print("Loading hotspot windows and WGS mutations...")
win_df = pd.read_csv(src_win)
snv_df = pd.read_csv(src_snv, low_memory=False)

snv_df["Chromosome"] = snv_df["Chromosome"].astype(str).str.replace("chr", "", case=False)
snv_df = snv_df[snv_df["Start_Position"].notna()].copy()
snv_df["Start_Position"] = snv_df["Start_Position"].astype(int)

n_pat = snv_df["patient"].nunique()

# 2. Local 100kb Regional Mutation Density Correction
print("Calculating local regional background mutation rates (100kb bins)...")
snv_df["bin_100k"] = (snv_df["Start_Position"] // 100000) * 100000
bin_counts = snv_df.groupby(["Chromosome", "bin_100k"]).size().to_dict()

chroms = win_df["Chromosome"].astype(str).str.replace("chr", "", case=False).to_numpy()
centers = win_df["Center_Position"].to_numpy(dtype=int)
n_obs_arr = win_df["Observed_Mutations"].to_numpy(dtype=int)

local_exps, local_obs_over_exps, nl10_ps_local = [], [], []

for c, pos, n_obs in zip(chroms, centers, n_obs_arr):
    bin_key = (c, (pos // 100000) * 100000)
    reg_muts = bin_counts.get(bin_key, 100)
    local_rate_per_bp = reg_muts / (100000.0 * n_pat)
    lam_local = max(n_pat * 50 * local_rate_per_bp, 0.05)

    local_exps.append(round(lam_local, 3))
    local_obs_over_exps.append(round(n_obs / max(lam_local, 1e-5), 2))
    nl10_ps_local.append(round(nb_nl10_sf(n_obs, 10.0, 10.0, lam_local), 3))

win_df["Local_Expected"] = local_exps
win_df["Local_Obs_Over_Exp"] = local_obs_over_exps
win_df["nl10_p_local"] = nl10_ps_local

win_df["nl10_q_local"] = bh_nlq(win_df["nl10_p_local"].to_numpy()).round(3)
win_df["q_value_local"] = np.power(10.0, -win_df["nl10_q_local"])
win_df = win_df.sort_values(["nl10_q_local", "Observed_Mutations"], ascending=False).reset_index(drop=True)

# 3. GENCODE Annotation with Index Mapping
gtf_path = find_gtf()
if gtf_path:
    print(f"Annotating windows using GENCODE GTF: {gtf_path.name}...")
    import gzip, re

    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    r_gn = re.compile(r'gene_name "([^"]+)"')
    r_gt = re.compile(r'gene_type "([^"]+)"|gene_biotype "([^"]+)"')

    gene_intervals = []
    with opener(gtf_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = f[0].replace("chr", "")
            start, end, strand = int(f[3]), int(f[4]), f[6]
            attr = f[8]
            gn = r_gn.search(attr)
            gt = r_gt.search(attr)
            gname = gn.group(1) if gn else "Unknown"
            gtype = (gt.group(1) or gt.group(2)) if gt else "Unknown"

            if strand == "+":
                p_start, p_end = max(1, start - 2000), start + 500
            else:
                p_start, p_end = max(1, end - 500), end + 2000

            gene_intervals.append({
                "chrom": chrom, "gname": gname, "gtype": gtype,
                "g_start": start, "g_end": end,
                "p_start": p_start, "p_end": p_end
            })

    genes_df = pd.DataFrame(gene_intervals)

    # Use index dictionary to guarantee row alignment
    win_df["chrom_clean"] = win_df["Chromosome"].astype(str).str.replace("chr", "", case=False)
    ann_dict = {}

    for c_clean, group_indices in win_df.groupby("chrom_clean", sort=False).groups.items():
        sub_win = win_df.loc[group_indices]
        sub_genes = genes_df[genes_df["chrom"] == c_clean]

        if sub_genes.empty:
            for idx in group_indices:
                ann_dict[idx] = ("Intergenic", "", "")
            continue

        p_starts = sub_genes["p_start"].to_numpy()
        p_ends = sub_genes["p_end"].to_numpy()
        g_starts = sub_genes["g_start"].to_numpy()
        g_ends = sub_genes["g_end"].to_numpy()
        gnames = sub_genes["gname"].to_numpy()
        gtypes = sub_genes["gtype"].to_numpy()

        for idx, pos in zip(group_indices, sub_win["Center_Position"]):
            p_hit = (pos >= p_starts) & (pos <= p_ends)
            if np.any(p_hit):
                ann_dict[idx] = ("Promoter", ";".join(np.unique(gnames[p_hit])[:3]), ";".join(np.unique(gtypes[p_hit])[:3]))
            else:
                b_hit = (pos >= g_starts) & (pos <= g_ends)
                if np.any(b_hit):
                    ann_dict[idx] = ("Gene Body / Intron", ";".join(np.unique(gnames[b_hit])[:3]), ";".join(np.unique(gtypes[b_hit])[:3]))
                else:
                    ann_dict[idx] = ("Intergenic", "", "")

    win_df["Annotation"] = [ann_dict[i][0] for i in win_df.index]
    win_df["Target_Gene"] = [ann_dict[i][1] for i in win_df.index]
    win_df["Gene_Type"] = [ann_dict[i][2] for i in win_df.index]
    win_df = win_df.drop(columns=["chrom_clean"])

# 4. Save Final Output
out_file = TAB / "25_prioritized_noncoding_hotspots.csv"
win_df.to_csv(out_file, index=False)
print(f"\nSaved prioritized noncoding hotspots to: {out_file}")

print("\n--- TOP 15 PRIORITIZED NONCODING HOTSPOTS ---")
display_cols = [c for c in ["Chromosome", "Center_Position", "Target_Gene", "Annotation", "Mutated_Patients", "Observed_Mutations", "Local_Expected", "q_value_local"] if c in win_df.columns]
print(win_df.head(15)[display_cols].to_string(index=False))