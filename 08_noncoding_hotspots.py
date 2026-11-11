"""
08_noncoding_hotspots.py

Fast, vectorized context-aware noncoding hotspot calling in PCAWG melanoma whole-genome data.

Inputs:
  - results/pcawg_melanoma_wgs_snvs.csv.gz (from script 07)

Outputs:
  - results/tables/23_noncoding_window_hotspots.csv
"""

import math
from pathlib import Path
import numpy as np
import pandas as pd

# ========================= CONFIGURATION ===========================
BASE = Path(r"D:\noncoding mutation finding surver")
RES = BASE / "results"
TAB = RES / "tables"
FIG = RES / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

WINDOW_SIZE = 50       # Window size in base pairs
MIN_PATIENTS = 3       # Minimum mutated donors to evaluate a window
FDR_THRESHOLD = 0.05
LN10 = math.log(10)
# ===================================================================


def nb_nl10_sf(k, shape, rate, lam0):
    """Negative binomial upper tail survival function in log10 space."""
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
    """Benjamini-Hochberg FDR adjustment (-log10 q-value)."""
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


# 1. Load & Sort WGS Mutations
src = RES / "pcawg_melanoma_wgs_snvs.csv.gz"
print(f"Loading melanoma WGS mutations from: {src.name}...")
df = pd.read_csv(src, low_memory=False)

df["Chromosome"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
df = df[df["Start_Position"].notna()].copy()
df["Start_Position"] = df["Start_Position"].astype(int)

# Sort mutations by chromosome and position for fast binary search
df = df.sort_values(["Chromosome", "Start_Position"]).reset_index(drop=True)

n_pat = df["patient"].nunique()
print(f"Loaded {len(df):,} mutations across {n_pat} melanoma patients.")

# 2. Global background rate estimation
genome_size_mb = 2800.0
mean_bg_rate = len(df) / (n_pat * genome_size_mb * 1e6)
print(f"Global noncoding mutation rate: {mean_bg_rate:.3e} mutations / bp / patient")

# 3. Candidate Recurrent Positions (>= 2 patients)
pos_counts = (
    df.groupby(["Chromosome", "Start_Position"])["patient"]
    .nunique()
    .reset_index(name="n_pts")
)
candidates = pos_counts[pos_counts["n_pts"] >= 2].copy()
print(f"Candidate recurrent positions (>= 2 patients): {len(candidates):,}")

# 4. Fast Window Evaluation via Binary Search
print("\nEvaluating 50bp candidate windows (fast binary search)...")
windows = []
seen_windows = set()

for chrom, sub_cands in candidates.groupby("Chromosome"):
    sub_all = df[df["Chromosome"] == chrom]
    all_pos = sub_all["Start_Position"].to_numpy()
    all_pts = sub_all["patient"].to_numpy()

    cand_positions = sub_cands["Start_Position"].to_numpy()

    for p in cand_positions:
        w_start = p - 25
        w_end = p + 25
        win_key = (chrom, w_start, w_end)

        if win_key in seen_windows:
            continue
        seen_windows.add(win_key)

        # Binary search range slice
        idx_start = np.searchsorted(all_pos, w_start, side="left")
        idx_end = np.searchsorted(all_pos, w_end, side="right")

        if idx_start >= idx_end:
            continue

        n_obs = idx_end - idx_start
        win_pts = len(np.unique(all_pts[idx_start:idx_end]))

        if win_pts >= MIN_PATIENTS:
            lam0 = n_pat * WINDOW_SIZE * mean_bg_rate
            alpha = 10.0
            nl10_p = nb_nl10_sf(n_obs, alpha, alpha, lam0)

            windows.append({
                "Chromosome": chrom,
                "Center_Position": p,
                "Window_Start": w_start,
                "Window_End": w_end,
                "Observed_Mutations": n_obs,
                "Mutated_Patients": win_pts,
                "Expected_Mutations": round(lam0, 3),
                "Obs_Over_Exp": round(n_obs / max(lam0, 1e-5), 2),
                "nl10_p": round(nl10_p, 3),
            })

win_df = pd.DataFrame(windows)
print(f"Total candidate windows tested: {len(win_df):,}")

if not win_df.empty:
    win_df["nl10_q"] = bh_nlq(win_df["nl10_p"].to_numpy()).round(3)
    win_df["q_value"] = np.power(10.0, -win_df["nl10_q"])
    win_df = win_df.sort_values(["nl10_q", "Observed_Mutations"], ascending=False)

    # Check TERT Promoter Positive Control
    tert_hits = win_df[
        (win_df["Chromosome"] == "5")
        & (win_df["Window_Start"].between(1295150, 1295300))
    ]

    print("\n--- SANITY CHECK: TERT Promoter Hotspot Detection ---")
    if not tert_hits.empty:
        print("TERT Promoter Window Results:")
        print(
            tert_hits[
                [
                    "Chromosome",
                    "Center_Position",
                    "Window_Start",
                    "Window_End",
                    "Mutated_Patients",
                    "Observed_Mutations",
                    "Expected_Mutations",
                    "q_value",
                ]
            ].to_string(index=False)
        )
    else:
        print("TERT promoter region did not meet the window filter thresholds.")

    # 5. Save Output
    out_win = TAB / "23_noncoding_window_hotspots.csv"
    win_df.to_csv(out_win, index=False)
    print(f"\nSaved noncoding hotspot window results to: {out_win}")

    sig_windows = win_df[win_df["q_value"] < FDR_THRESHOLD]
    print(
        f"Statistically significant noncoding hotspot windows (q < {FDR_THRESHOLD}): {len(sig_windows):,}"
    )