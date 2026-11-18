# Experiment 6: gnomAD allele-frequency spectrum of NC-MutScan v3.3 calls

- Build hash recorded: `not given`; engine `my_engine:score`; flank ±25 bp; seed 42.
- gnomAD files: gnomad.genomes.v4.1.sites.chr21.vcf.bgz; common cut-off AF ≥ 0.01; PASS biallelic SNVs, AN ≥ 80% of the 90th-percentile AN; CDS of any isoform excluded.
- Domains: utr (exon, non-CDS), splice (intronic 1–100 bp from an exon edge), promoter (TSS −1000/+100). Scored sample: 92,580 variants.
- Strata: domain-distance bin × CpG/ts-tv (4) × GC bin (4). Intervals: block bootstrap, 100 kb blocks, 200 replicates (50 for the slope).

## Variants per frequency class

| Class | Definition | Seen in domains | Sampled |
|---|---|---|---|
| singleton | AC = 1 | 46919 | 46919 |
| ultra_rare | AC >= 2, AF < 1e-4 | 34145 | 34145 |
| rare | 1e-4 <= AF < 1e-3 | 6472 | 6472 |
| low_freq | 1e-3 <= AF < common cut | 2891 | 2891 |
| common | common cut <= AF < 5% | 1049 | 1049 |
| very_common | AF >= 5% | 1104 | 1104 |

*O/E = observed calls / calls expected from the stratum-specific call rate (1 = no frequency effect). MH OR = Mantel–Haenszel odds ratio of a call versus the **common** class. Slope = stratified log-odds change of a call per frequency-class step (negative = calls fall as frequency rises).*

## call | all

Trend: slope +0.010 per class step (95% CI -0.022 to +0.050; bootstrap two-sided P = 0.400).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 46,919 | 2,394 | 2398.7 | 1.00 (0.97–1.02) | 0.94 (0.75–1.31) |
| ultra_rare | 34,145 | 1,664 | 1665.5 | 1.00 (0.96–1.03) | 0.96 (0.76–1.28) |
| rare | 6,472 | 298 | 303.8 | 0.98 (0.86–1.08) | 0.93 (0.71–1.31) |
| low_freq | 2,891 | 140 | 136.8 | 1.02 (0.91–1.17) | 0.98 (0.76–1.29) |
| common | 1,049 | 53 | 50.6 | 1.05 (0.79–1.28) | reference |
| very_common | 1,104 | 59 | 52.7 | 1.12 (0.86–1.42) | 1.01 (0.65–1.59) |

## call | splice

Trend: slope -0.022 per class step (95% CI -0.080 to +0.049; bootstrap two-sided P = 0.360).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 10,993 | 660 | 653.4 | 1.01 (0.94–1.05) | 1.45 (0.73–3.52) |
| ultra_rare | 7,892 | 443 | 442.7 | 1.00 (0.94–1.06) | 1.43 (0.81–3.43) |
| rare | 1,492 | 74 | 77.4 | 0.96 (0.79–1.17) | 1.23 (0.64–3.52) |
| low_freq | 679 | 34 | 34.2 | 0.99 (0.70–1.29) | 1.50 (0.61–4.67) |
| common | 241 | 10 | 13.4 | 0.75 (0.34–1.22) | reference |
| very_common | 285 | 15 | 14.9 | 1.01 (0.65–1.50) | 1.33 (0.54–4.69) |

## call | utr

Trend: slope +0.013 per class step (95% CI -0.039 to +0.055; bootstrap two-sided P = 0.760).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 20,679 | 1,018 | 1004.7 | 1.01 (0.98–1.05) | 0.86 (0.62–1.69) |
| ultra_rare | 14,880 | 679 | 699.0 | 0.97 (0.93–1.01) | 0.83 (0.58–1.54) |
| rare | 2,602 | 113 | 119.2 | 0.95 (0.76–1.15) | 0.81 (0.45–1.66) |
| low_freq | 1,128 | 56 | 52.2 | 1.07 (0.81–1.34) | 0.90 (0.62–1.60) |
| common | 407 | 22 | 18.9 | 1.17 (0.63–1.62) | reference |
| very_common | 476 | 28 | 22.2 | 1.26 (0.86–1.70) | 1.05 (0.54–2.44) |

## call | promoter

Trend: slope +0.029 per class step (95% CI -0.005 to +0.081; bootstrap two-sided P = 0.160).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 15,247 | 716 | 740.7 | 0.97 (0.91–1.01) | 0.82 (0.62–1.21) |
| ultra_rare | 11,373 | 542 | 523.8 | 1.03 (0.97–1.10) | 0.91 (0.71–1.34) |
| rare | 2,378 | 111 | 107.2 | 1.04 (0.90–1.23) | 0.94 (0.68–1.33) |
| low_freq | 1,084 | 50 | 50.4 | 0.99 (0.80–1.30) | 0.86 (0.51–1.52) |
| common | 401 | 21 | 18.3 | 1.15 (0.79–1.47) | reference |
| very_common | 343 | 16 | 15.7 | 1.02 (0.60–1.49) | 0.80 (0.42–1.68) |

## call | canonical_1-2

Trend: slope +0.076 per class step (95% CI -0.132 to +0.248; bootstrap two-sided P = 0.480).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 256 | 119 | 121.7 | 0.98 (0.91–1.06) | 0.65 (0.20–3.81) |
| ultra_rare | 167 | 78 | 77.8 | 1.00 (0.90–1.11) | 0.69 (0.16–3.89) |
| rare | 22 | 12 | 9.7 | 1.23 (0.63–1.58) | 0.69 (0.19–4.81) |
| low_freq | 9 | 3 | 3.4 | 0.89 (0.21–1.90) | NA (NA–NA) |
| common | 6 | 3 | 2.6 | 1.17 (0.00–1.99) | reference |
| very_common | 4 | 2 | 1.9 | 1.07 (0.00–2.17) | 2.00 (0.60–3.61) |

## call | splice d3-5

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 298 | 13 | 14.4 | 0.91 (0.62–1.22) | NA (NA–NA) |
| ultra_rare | 233 | 11 | 9.6 | 1.14 (0.69–1.61) | NA (NA–NA) |
| rare | 32 | 1 | 1.1 | 0.94 (0.00–2.79) | NA (NA–NA) |
| low_freq | 15 | 0 | 0.6 | 0.00 (0.00–0.00) | NA (NA–NA) |
| common | 3 | 0 | 0.1 | 0.00 (0.00–0.00) | reference |
| very_common | 6 | 1 | 0.3 | 3.95 (0.00–8.22) | NA (NA–NA) |

## call | splice d6-10

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 530 | 33 | 28.0 | 1.18 (0.95–1.41) | NA (NA–NA) |
| ultra_rare | 364 | 16 | 17.6 | 0.91 (0.58–1.24) | NA (NA–NA) |
| rare | 74 | 2 | 3.8 | 0.53 (0.00–1.54) | NA (NA–NA) |
| low_freq | 43 | 0 | 2.2 | 0.00 (0.00–0.00) | NA (NA–NA) |
| common | 9 | 0 | 0.4 | 0.00 (0.00–0.00) | reference |
| very_common | 20 | 2 | 0.9 | 2.14 (0.00–5.63) | NA (NA–NA) |

## call | splice d11-20

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 1,125 | 53 | 57.5 | 0.92 (0.74–1.10) | NA (NA–NA) |
| ultra_rare | 827 | 45 | 40.8 | 1.10 (0.89–1.34) | NA (NA–NA) |
| rare | 152 | 8 | 7.0 | 1.14 (0.37–1.86) | NA (NA–NA) |
| low_freq | 67 | 4 | 3.3 | 1.22 (0.36–2.23) | NA (NA–NA) |
| common | 22 | 0 | 1.2 | 0.00 (0.00–0.00) | reference |
| very_common | 29 | 1 | 1.3 | 0.78 (0.00–2.49) | NA (NA–NA) |

## call | splice d21-50

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 3,442 | 161 | 157.8 | 1.02 (0.91–1.11) | 1.79 (0.61–4.44) |
| ultra_rare | 2,431 | 108 | 106.7 | 1.01 (0.88–1.18) | 1.71 (0.66–4.40) |
| rare | 491 | 22 | 21.4 | 1.03 (0.74–1.40) | 1.61 (0.59–4.70) |
| low_freq | 213 | 5 | 9.0 | 0.56 (0.00–1.06) | 0.67 (0.07–2.11) |
| common | 78 | 2 | 3.4 | 0.59 (0.00–1.56) | reference |
| very_common | 84 | 4 | 3.7 | 1.09 (0.27–2.27) | 1.31 (0.21–7.89) |

## call | splice d51-100

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 5,342 | 281 | 274.0 | 1.03 (0.94–1.10) | 1.24 (0.55–5.11) |
| ultra_rare | 3,870 | 185 | 190.1 | 0.97 (0.87–1.09) | 1.14 (0.53–4.64) |
| rare | 721 | 29 | 34.3 | 0.84 (0.58–1.09) | 0.95 (0.37–3.81) |
| low_freq | 332 | 22 | 15.8 | 1.40 (1.00–1.96) | 1.65 (0.70–7.57) |
| common | 123 | 5 | 5.8 | 0.86 (0.15–1.74) | reference |
| very_common | 142 | 5 | 6.9 | 0.73 (0.17–1.44) | 0.87 (0.25–3.32) |

## s_raw>=2.8 | all

Trend: slope +0.008 per class step (95% CI -0.019 to +0.044; bootstrap two-sided P = 0.560).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 46,919 | 4,678 | 4668.8 | 1.00 (0.99–1.02) | 0.99 (0.83–1.19) |
| ultra_rare | 34,145 | 3,279 | 3313.0 | 0.99 (0.97–1.01) | 0.97 (0.82–1.19) |
| rare | 6,472 | 616 | 610.9 | 1.01 (0.95–1.07) | 0.99 (0.84–1.22) |
| low_freq | 2,891 | 281 | 274.5 | 1.02 (0.91–1.14) | 1.01 (0.82–1.29) |
| common | 1,049 | 102 | 100.4 | 1.02 (0.86–1.18) | reference |
| very_common | 1,104 | 118 | 106.4 | 1.11 (0.90–1.35) | 1.12 (0.82–1.49) |

## s_raw>=4.2 | all

Trend: slope +0.003 per class step (95% CI -0.037 to +0.033; bootstrap two-sided P = 0.800).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 46,919 | 1,796 | 1772.4 | 1.01 (0.98–1.05) | 0.92 (0.70–1.36) |
| ultra_rare | 34,145 | 1,210 | 1238.7 | 0.98 (0.94–1.01) | 0.90 (0.69–1.34) |
| rare | 6,472 | 222 | 226.2 | 0.98 (0.83–1.10) | 0.88 (0.63–1.28) |
| low_freq | 2,891 | 100 | 101.0 | 0.99 (0.81–1.17) | 0.91 (0.69–1.28) |
| common | 1,049 | 41 | 37.7 | 1.09 (0.75–1.40) | reference |
| very_common | 1,104 | 46 | 39.0 | 1.18 (0.93–1.62) | 1.04 (0.68–1.82) |

## s_raw>=6 | all

Trend: slope +0.009 per class step (95% CI -0.053 to +0.063; bootstrap two-sided P = 0.560).

| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |
|---|---|---|---|---|---|
| singleton | 46,919 | 1,012 | 991.6 | 1.02 (0.98–1.06) | 0.92 (0.66–1.58) |
| ultra_rare | 34,145 | 676 | 706.7 | 0.96 (0.90–1.02) | 0.89 (0.61–1.61) |
| rare | 6,472 | 138 | 130.7 | 1.06 (0.91–1.21) | 0.94 (0.65–1.65) |
| low_freq | 2,891 | 50 | 58.6 | 0.85 (0.60–1.08) | 0.83 (0.53–1.56) |
| common | 1,049 | 24 | 22.1 | 1.08 (0.65–1.49) | reference |
| very_common | 1,104 | 32 | 22.3 | 1.44 (1.01–1.81) | 1.24 (0.76–2.97) |

## How to read this

- Classes are sampled at different rates, so O/E is relative to the sampled class mix: use its pattern across classes, and take effect sizes from the MH OR and the slope.
- Expected pattern if calls mark functional disruption: O/E > 1 for singletons, falling below 1 for the common classes, strongest for canonical ±1–2 bp splice sites and at the higher score thresholds (dose-response).
- No frequency gradient is itself informative: it limits the claim that calls mark selected positions in that domain (for example promoter calls from generic motifs).
- Remaining confounders: CpG methylation level, local recombination/background selection and gene-level constraint are not matched; report them as limitations.