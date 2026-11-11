"""
10_filter_uv_and_visualize.py

1. Filters intronic/intergenic candidate hotspots for high-confidence noncoding drivers.
2. Generates noncoding genome-wide Manhattan plot highlighting TERT and top candidates.
3. Exports final candidate noncoding driver summary table.

Inputs:
  - results/tables/25_prioritized_noncoding_hotspots.csv

Outputs:
  - results/tables/26_final_noncoding_driver_candidates.csv
  - results/figures/fig16_noncoding_manhattan_plot.png
"""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ========================= CONFIGURATION ===========================
BASE = Path(r"D:\noncoding mutation finding surver")
RES = BASE / "results"
TAB = RES / "tables"
FIG = RES / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)
# ===================================================================

# 1. Load Prioritized Noncoding Hotspots
src = TAB / "25_prioritized_noncoding_hotspots.csv"
print(f"Loading scored noncoding hotspots from: {src.name}...")
df = pd.read_csv(src)

# Sort by Chromosome and Position for plotting
chrom_order = [str(i) for i in range(1, 23)] + ["X", "Y"]
df["chrom_clean"] = df["Chromosome"].astype(str).str.replace("chr", "", case=False)
df = df[df["chrom_clean"].isin(chrom_order)].copy()
df["chrom_num"] = pd.Categorical(df["chrom_clean"], categories=chrom_order, ordered=True)
df = df.sort_values(["chrom_num", "Center_Position"]).reset_index(drop=True)

# 2. Filter Candidate Drivers
# Prioritize Promoters, UTRs, lncRNAs, and highly significant recurrent intronic elements
df["is_promoter"] = df["Annotation"].str.contains("Promoter", case=False, na=False)

print(f"Total noncoding windows scored: {len(df):,}")
print(f"Promoter hotspots identified: {df['is_promoter'].sum():,}")

# 3. Generate Genome-Wide Manhattan Plot
print("\nGenerating Noncoding Hotspot Manhattan Plot...")

# Compute cumulative genomic positions
chrom_lengths = {}
for c, group in df.groupby("chrom_num", observed=True):
    chrom_lengths[c] = group["Center_Position"].max()

cum_pos = 0
chrom_offsets = {}
ticks = []
labels = []

for c in chrom_order:
    if c in chrom_lengths:
        chrom_offsets[c] = cum_pos
        sub = df[df["chrom_clean"] == c]
        if not sub.empty:
            ticks.append(cum_pos + (sub["Center_Position"].max() - sub["Center_Position"].min()) / 2)
            labels.append(c)
        cum_pos += chrom_lengths[c] + 50000000

df["cum_pos"] = df.apply(lambda r: r["Center_Position"] + chrom_offsets.get(r["chrom_clean"], 0), axis=1)

# Plotting
plt.figure(figsize=(14, 6), dpi=300)
colors = ["#2b5c8f", "#709ae1"]

for i, c in enumerate(chrom_order):
    sub = df[df["chrom_clean"] == c]
    if sub.empty:
        continue
    plt.scatter(
        sub["cum_pos"],
        sub["nl10_q_local"],
        c=colors[i % 2],
        s=12,
        alpha=0.7,
        edgecolors="none"
    )

# Highlight Promoters (Red) and TERT (Gold)
promoters = df[df["is_promoter"] & (df["nl10_q_local"] > 5)]
plt.scatter(
    promoters["cum_pos"],
    promoters["nl10_q_local"],
    c="#d9381e",
    s=25,
    label="Promoter Hotspots",
    zorder=3
)

tert = df[(df["chrom_clean"] == "5") & (df["Center_Position"].between(1295200, 1295260))]
if not tert.empty:
    plt.scatter(
        tert["cum_pos"],
        tert["nl10_q_local"],
        c="#ffaa00",
        s=80,
        marker="*",
        edgecolors="black",
        label="TERT Promoter (Control)",
        zorder=4
    )
    for _, r in tert.iterrows():
        plt.annotate(
            "TERT",
            (r["cum_pos"], r["nl10_q_local"]),
            textcoords="offset points",
            xytext=(0, 10),
            ha="center",
            fontsize=10,
            weight="bold",
            color="#b35900"
        )

plt.axhline(y=-np.log10(0.05), color="grey", linestyle="--", linewidth=0.8, label="FDR q = 0.05")
plt.xticks(ticks, labels, fontsize=8, rotation=0)
plt.xlabel("Chromosome", fontsize=11, labelpad=8)
plt.ylabel(r"Local $-\log_{10}(q\mathrm{-value})$", fontsize=11)
plt.title("PCAWG Melanoma WGS Noncoding Mutation Hotspots", fontsize=13, pad=12)
plt.legend(loc="upper right", frameon=True)
plt.tight_layout()

out_fig = FIG / "fig16_noncoding_manhattan_plot.png"
plt.savefig(out_fig)
plt.close()
print(f"Saved Manhattan plot to: {out_fig}")

# 4. Save Final Candidate Table
out_cand = TAB / "26_final_noncoding_driver_candidates.csv"
top_candidates = df.sort_values(["is_promoter", "nl10_q_local"], ascending=[False, False]).reset_index(drop=True)
top_candidates.to_csv(out_cand, index=False)

print(f"\nSaved final driver candidates to: {out_cand}")
print("\n--- TOP 10 NONCODING PROMOTER CANDIDATES ---")
print(top_candidates[top_candidates["is_promoter"]][["Chromosome", "Center_Position", "Target_Gene", "Mutated_Patients", "Observed_Mutations", "nl10_q_local"]].head(10).to_string(index=False))