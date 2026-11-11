import csv
import os
from collections import defaultdict, Counter

INPUT_FILE = "./results/pan_cancer_noncoding_mutations.csv"
OUTPUT_HOTSPOTS = "./results/noncoding_hotspots_summary.csv"
WINDOW_SIZE = 50

print("[*] Starting non-coding hotspot analysis...", flush=True)

if not os.path.exists(INPUT_FILE):
    print(f"[!] ERROR: Input file '{INPUT_FILE}' not found! Verify step 16 generated it.", flush=True)
    exit(1)

hotspots = defaultdict(list)
gene_counts = Counter()
class_counts = Counter()
total_rows = 0

with open(INPUT_FILE, 'r', encoding='utf-8', errors='ignore') as f:
    reader = csv.DictReader(f)
    for row in reader:
        total_rows += 1
        chrom = row.get('Chromosome', '')
        try:
            start = int(row.get('Start_Position', 0))
        except (ValueError, TypeError):
            continue
        
        gene = row.get('Hugo_Symbol', 'Intergenic/Unknown')
        v_class = row.get('Variant_Classification', 'Unknown')
        cancer_type = row.get('Cancer_Type', 'Unknown')
        sample = row.get('Tumor_Sample_Barcode', '')

        # Bin mutations into 50bp genomic windows
        window_bin = (start // WINDOW_SIZE) * WINDOW_SIZE
        key = (chrom, window_bin, gene)

        hotspots[key].append({
            'start': start,
            'v_class': v_class,
            'cancer_type': cancer_type,
            'sample': sample
        })
        
        if gene:
            gene_counts[gene] += 1
        class_counts[v_class] += 1

print(f"[+] Loaded {total_rows} variants across {len(gene_counts)} gene/regulatory regions.", flush=True)

OUTPUT_FIELDS = [
    'Chromosome', 'Window_Start', 'Window_End', 'Gene', 
    'Recurrence_Count', 'Sample_Count', 'Cancer_Types', 'Variant_Types'
]

recurrent_hotspots = 0

with open(OUTPUT_HOTSPOTS, 'w', newline='', encoding='utf-8') as out_f:
    writer = csv.DictWriter(out_f, fieldnames=OUTPUT_FIELDS)
    writer.writeheader()

    sorted_keys = sorted(hotspots.keys(), key=lambda k: len(hotspots[k]), reverse=True)

    for key in sorted_keys:
        mutations = hotspots[key]
        rec_count = len(mutations)
        
        if rec_count >= 2:
            recurrent_hotspots += 1

        chrom, win_start, gene = key
        cancers = ";".join(sorted(set(m['cancer_type'] for m in mutations if m['cancer_type'])))
        unique_samples = len(set(m['sample'] for m in mutations if m['sample']))
        v_types = ";".join(sorted(set(m['v_class'] for m in mutations if m['v_class'])))

        writer.writerow({
            'Chromosome': chrom,
            'Window_Start': win_start,
            'Window_End': win_start + WINDOW_SIZE,
            'Gene': gene,
            'Recurrence_Count': rec_count,
            'Sample_Count': unique_samples,
            'Cancer_Types': cancers,
            'Variant_Types': v_types
        })

print(f"\n[+] SUCCESS: Saved hotspot analysis to: {OUTPUT_HOTSPOTS}", flush=True)
print(f"[+] Recurrent non-coding hotspot regions found: {recurrent_hotspots}", flush=True)

print("\n--- Non-Coding Variant Types Breakdown ---", flush=True)
for vtype, count in class_counts.most_common(8):
    print(f"  - {vtype}: {count}", flush=True)

print("\n--- Top 10 Targeted Non-Coding Genes/Loci ---", flush=True)
for gene, count in gene_counts.most_common(10):
    if gene and gene != 'Intergenic/Unknown':
        print(f"  - {gene}: {count} variants", flush=True)