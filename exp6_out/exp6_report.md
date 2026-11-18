# Experiment 6: gnomAD allele-frequency spectrum of NC-MutScan v3.3 calls

- Build hash recorded: `b276d0d44045`; engine `my_engine:score`; flank ±25 bp; seed 42.
- gnomAD files: gnomad.genomes.v4.1.sites.chr21.vcf.bgz, gnomad.genomes.v4.1.sites.chr22.vcf.bgz; common cut-off AF ≥ 0.01; PASS biallelic SNVs, AN ≥ 80% of the 90th-percentile AN; CDS of any isoform excluded.
- Domains: utr (exon, non-CDS), splice (intronic 1–100 bp from an exon edge), promoter (TSS −1000/+100). Scored sample: 503,112 variants.
- Strata: domain-distance bin × CpG/ts-tv (4) × GC bin (4). Intervals: block bootstrap, 100 kb blocks, 1000 replicates (200 for the slope).

## Variants per frequency class

| Class | Definition | Seen in domains | Sampled |
|---|---|---|---|
| singleton | AC = 1 | 897078 | 150000 |
| ultra_rare | AC >= 2, AF < 1e-4 | 671089 | 150000 |
| rare | 1e-4 <= AF < 1e-3 | 117097 | 117097 |
| low_freq | 1e-3 <= AF < common cut | 45603 | 45603 |
| common | common cut <= AF < 5% | 17889 | 17889 |
| very_common | AF >= 5% | 22523 | 22523 |

*O/E = observed calls / calls expected from the stratum-specific call rate (1 = no frequency effect). MH OR = Mantel–Haenszel odds ratio of a call versus the **common** class. Slope = stratified log-odds change of a call per frequency-class step (negative = calls fall as frequency rises).*

## call | all

Trend: slope -0.002 per class step (95% CI -0.011 to +0.008; bootstrap two-sided P = 0.730).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 150,000 | 7,256 | 7210.7 | 1.01 (0.99–1.03) | 0.96 (0.89–1.04) |
| ultra_rare | 150,000 | 6,834 | 6814.6 | 1.00 (0.98–1.02) | 0.97 (0.90–1.04) |
| rare | 117,097 | 5,121 | 5172.0 | 0.99 (0.97–1.01) | 0.96 (0.88–1.03) |
| low_freq | 45,603 | 1,971 | 2029.0 | 0.97 (0.93–1.01) | 0.94 (0.87–1.02) |
| common | 17,889 | 825 | 797.7 | 1.03 (0.97–1.10) | reference |
| very_common | 22,523 | 1,027 | 1010.0 | 1.02 (0.95–1.08) | 0.98 (0.89–1.09) |

## call | splice

Trend: slope -0.000 per class step (95% CI -0.019 to +0.023; bootstrap two-sided P = 0.840).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 34,843 | 1,867 | 1811.6 | 1.03 (1.00–1.06) | 1.04 (0.89–1.24) |
| ultra_rare | 34,728 | 1,615 | 1673.0 | 0.97 (0.93–1.00) | 1.01 (0.86–1.20) |
| rare | 26,762 | 1,248 | 1229.8 | 1.01 (0.97–1.06) | 1.05 (0.90–1.26) |
| low_freq | 10,375 | 433 | 476.9 | 0.91 (0.83–0.99) | 0.95 (0.78–1.17) |
| common | 4,045 | 179 | 186.1 | 0.96 (0.83–1.10) | reference |
| very_common | 4,994 | 265 | 229.5 | 1.15 (1.04–1.28) | 1.22 (0.99–1.52) |

## call | utr

Trend: slope -0.006 per class step (95% CI -0.021 to +0.008; bootstrap two-sided P = 0.490).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 69,163 | 3,283 | 3263.4 | 1.01 (0.98–1.03) | 0.98 (0.87–1.10) |
| ultra_rare | 70,021 | 3,154 | 3118.2 | 1.01 (0.98–1.04) | 0.99 (0.88–1.11) |
| rare | 54,560 | 2,326 | 2361.9 | 0.98 (0.95–1.02) | 0.96 (0.86–1.08) |
| low_freq | 21,103 | 906 | 923.6 | 0.98 (0.92–1.04) | 0.96 (0.85–1.10) |
| common | 8,250 | 371 | 362.9 | 1.02 (0.92–1.13) | reference |
| very_common | 10,246 | 445 | 455.0 | 0.98 (0.89–1.06) | 0.95 (0.82–1.10) |

## call | promoter

Trend: slope +0.004 per class step (95% CI -0.014 to +0.020; bootstrap two-sided P = 0.830).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 45,994 | 2,106 | 2135.7 | 0.99 (0.95–1.02) | 0.89 (0.78–1.04) |
| ultra_rare | 45,251 | 2,065 | 2023.3 | 1.02 (0.99–1.06) | 0.92 (0.82–1.07) |
| rare | 35,775 | 1,547 | 1580.2 | 0.98 (0.94–1.02) | 0.88 (0.78–1.02) |
| low_freq | 14,125 | 632 | 628.4 | 1.01 (0.93–1.08) | 0.91 (0.79–1.06) |
| common | 5,594 | 275 | 248.7 | 1.11 (0.97–1.23) | reference |
| very_common | 7,283 | 317 | 325.5 | 0.97 (0.87–1.08) | 0.89 (0.74–1.06) |

## call | canonical_1-2

Trend: slope -0.084 per class step (95% CI -0.167 to +0.016; bootstrap two-sided P = 0.100).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 637 | 280 | 256.4 | 1.09 (1.01–1.17) | 1.55 (0.78–3.25) |
| ultra_rare | 550 | 189 | 205.5 | 0.92 (0.83–1.01) | 1.20 (0.60–2.44) |
| rare | 355 | 125 | 125.6 | 0.99 (0.87–1.11) | 1.28 (0.64–2.68) |
| low_freq | 127 | 44 | 43.6 | 1.01 (0.78–1.27) | 1.64 (0.73–3.70) |
| common | 50 | 14 | 17.0 | 0.83 (0.50–1.22) | reference |
| very_common | 49 | 15 | 18.9 | 0.79 (0.50–1.13) | 0.81 (0.23–2.05) |

## call | splice d3-5

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 1,084 | 42 | 40.3 | 1.04 (0.82–1.28) | 0.53 (0.24–1.78) |
| ultra_rare | 1,031 | 29 | 34.3 | 0.84 (0.58–1.11) | 0.44 (0.19–1.48) |
| rare | 706 | 23 | 23.7 | 0.97 (0.64–1.31) | 0.53 (0.23–1.77) |
| low_freq | 243 | 9 | 7.6 | 1.18 (0.53–1.93) | 0.73 (0.21–2.64) |
| common | 115 | 7 | 3.9 | 1.79 (0.60–3.11) | reference |
| very_common | 124 | 4 | 4.1 | 0.98 (0.21–2.09) | 0.60 (0.09–2.18) |

## call | splice d6-10

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 1,737 | 65 | 62.8 | 1.03 (0.85–1.21) | 1.38 (0.67–4.89) |
| ultra_rare | 1,837 | 59 | 61.0 | 0.97 (0.76–1.18) | 1.23 (0.58–4.47) |
| rare | 1,353 | 45 | 44.1 | 1.02 (0.75–1.28) | 1.26 (0.59–4.20) |
| low_freq | 514 | 13 | 16.6 | 0.78 (0.42–1.22) | 1.04 (0.40–4.01) |
| common | 229 | 6 | 7.5 | 0.79 (0.24–1.49) | reference |
| very_common | 236 | 12 | 7.9 | 1.52 (0.68–2.40) | 2.12 (0.72–8.65) |

## call | splice d11-20

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 3,669 | 162 | 147.0 | 1.10 (0.98–1.24) | 0.97 (0.61–1.76) |
| ultra_rare | 3,725 | 139 | 141.3 | 0.98 (0.84–1.13) | 0.80 (0.49–1.47) |
| rare | 2,962 | 98 | 108.4 | 0.90 (0.75–1.06) | 0.75 (0.46–1.35) |
| low_freq | 1,168 | 38 | 43.6 | 0.87 (0.62–1.14) | 0.70 (0.40–1.39) |
| common | 391 | 17 | 14.4 | 1.18 (0.69–1.76) | reference |
| very_common | 570 | 22 | 21.3 | 1.03 (0.63–1.50) | 0.88 (0.47–1.70) |

## call | splice d21-50

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 10,932 | 504 | 500.0 | 1.01 (0.94–1.08) | 0.99 (0.75–1.34) |
| ultra_rare | 10,885 | 440 | 475.0 | 0.93 (0.85–1.00) | 0.94 (0.71–1.29) |
| rare | 8,442 | 407 | 358.2 | 1.14 (1.05–1.23) | 1.16 (0.86–1.62) |
| low_freq | 3,352 | 112 | 143.3 | 0.78 (0.65–0.93) | 0.79 (0.57–1.13) |
| common | 1,231 | 52 | 52.9 | 0.98 (0.74–1.26) | reference |
| very_common | 1,616 | 84 | 69.6 | 1.21 (0.98–1.42) | 1.24 (0.88–1.71) |

## call | splice d51-100

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 16,784 | 814 | 805.0 | 1.01 (0.95–1.07) | 1.06 (0.84–1.38) |
| ultra_rare | 16,700 | 759 | 755.9 | 1.00 (0.94–1.06) | 1.10 (0.89–1.40) |
| rare | 12,944 | 550 | 569.7 | 0.97 (0.89–1.04) | 1.05 (0.84–1.37) |
| low_freq | 4,971 | 217 | 222.3 | 0.98 (0.85–1.10) | 1.06 (0.83–1.42) |
| common | 2,029 | 83 | 90.4 | 0.92 (0.73–1.11) | reference |
| very_common | 2,399 | 128 | 107.7 | 1.19 (1.00–1.37) | 1.31 (0.99–1.74) |

## s_raw>=2.8 | all

Trend: slope -0.001 per class step (95% CI -0.007 to +0.006; bootstrap two-sided P = 0.770).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 150,000 | 14,562 | 14508.5 | 1.00 (0.99–1.02) | 0.97 (0.92–1.02) |
| ultra_rare | 150,000 | 14,026 | 14003.4 | 1.00 (0.99–1.01) | 0.97 (0.92–1.02) |
| rare | 117,097 | 10,635 | 10724.6 | 0.99 (0.98–1.01) | 0.96 (0.91–1.01) |
| low_freq | 45,603 | 4,167 | 4198.5 | 0.99 (0.97–1.02) | 0.96 (0.91–1.02) |
| common | 17,889 | 1,700 | 1651.8 | 1.03 (0.98–1.07) | reference |
| very_common | 22,523 | 2,092 | 2095.3 | 1.00 (0.96–1.04) | 0.97 (0.91–1.04) |

## s_raw>=4.2 | all

Trend: slope +0.002 per class step (95% CI -0.010 to +0.014; bootstrap two-sided P = 0.830).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 150,000 | 5,265 | 5245.1 | 1.00 (0.98–1.03) | 0.96 (0.88–1.05) |
| ultra_rare | 150,000 | 5,077 | 5075.5 | 1.00 (0.98–1.02) | 0.96 (0.89–1.06) |
| rare | 117,097 | 3,877 | 3895.3 | 1.00 (0.97–1.02) | 0.96 (0.88–1.05) |
| low_freq | 45,603 | 1,460 | 1520.7 | 0.96 (0.91–1.01) | 0.93 (0.84–1.03) |
| common | 17,889 | 619 | 595.7 | 1.04 (0.96–1.12) | reference |
| very_common | 22,523 | 783 | 748.8 | 1.05 (0.97–1.12) | 1.01 (0.89–1.12) |

## s_raw>=6 | all

Trend: slope -0.002 per class step (95% CI -0.019 to +0.013; bootstrap two-sided P = 0.870).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 150,000 | 3,330 | 3280.0 | 1.02 (0.99–1.04) | 0.98 (0.88–1.11) |
| ultra_rare | 150,000 | 3,165 | 3203.8 | 0.99 (0.96–1.01) | 0.96 (0.87–1.07) |
| rare | 117,097 | 2,471 | 2460.4 | 1.00 (0.97–1.04) | 0.98 (0.88–1.09) |
| low_freq | 45,603 | 910 | 957.3 | 0.95 (0.89–1.01) | 0.93 (0.82–1.05) |
| common | 17,889 | 385 | 374.4 | 1.03 (0.93–1.13) | reference |
| very_common | 22,523 | 485 | 470.1 | 1.03 (0.94–1.12) | 1.01 (0.87–1.17) |

## How to read this

- Classes are sampled at different rates, so O/E is relative to the sampled class mix: use its pattern across classes, and take effect sizes from the MH OR and the slope.
- Expected pattern if calls mark functional disruption: O/E > 1 for singletons, falling below 1 for the common classes, strongest for canonical ±1–2 bp splice sites and at the higher score thresholds (dose-response).
- No frequency gradient is itself informative: it limits the claim that calls mark selected positions in that domain (for example promoter calls from generic motifs).
- Remaining confounders: CpG methylation level, local recombination/background selection and gene-level constraint are not matched; report them as limitations.