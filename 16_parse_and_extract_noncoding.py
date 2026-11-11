import os
import csv
import gzip

DATA_DIR = "./maf_files"
OUTPUT_FILE = "./results/pan_cancer_noncoding_mutations.csv"

# Non-coding region classifications
NONCODING_KEYWORDS = {
    "5'flank", "3'flank", "5'utr", "3'utr", "intron", "igr", 
    "rna", "splice_site", "lincrna", "targeted_region", "promoter",
    "upstream", "downstream", "intergenic", "splice", "non_coding"
}

SKIP_FILES = {'annotations.txt', 'manifest.txt', 'sample_sheet.tsv'}

OUTPUT_FIELDS = [
    'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 
    'Tumor_Seq_Allele2', 'Variant_Classification', 'Hugo_Symbol', 
    'Tumor_Sample_Barcode', 'Cancer_Type'
]

os.makedirs("./results", exist_ok=True)

print("[*] Scanning `./maf_files` for mutation callsets (including compressed .gz files)...")

mutation_files = []
for root, dirs, files in os.walk(DATA_DIR):
    for f in files:
        f_lower = f.lower()
        if f_lower in SKIP_FILES:
            continue
        # Skip logs, metadata, and parcel files
        if any(skip in f_lower for skip in ['clinical', 'meta_', 'case_list', 'patient', 'sample_info', '.parcel', 'logs']):
            continue
        
        # Capture MAF, MAF.GZ, TSV, TSV.GZ, and DONOR files
        if f_lower.endswith(('.maf', '.maf.gz', '.tsv', '.tsv.gz', '.txt', '.txt.gz', '.donor')):
            mutation_files.append(os.path.join(root, f))

print(f"[*] Identified {len(mutation_files)} mutation callset file(s).")

total_count = 0

with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8') as out_f:
    writer = csv.DictWriter(out_f, fieldnames=OUTPUT_FIELDS, extrasaction='ignore')
    writer.writeheader()

    for file in mutation_files:
        cancer_type = os.path.basename(os.path.dirname(file))
        print(f"\n[>] Extracting variants from: {file}")
        
        try:
            # Handle gzipped (.gz) vs uncompressed text files dynamically
            if file.endswith('.gz'):
                f_in = gzip.open(file, 'rt', encoding='utf-8', errors='ignore')
            else:
                f_in = open(file, 'r', encoding='utf-8', errors='ignore')

            with f_in as f:
                # Filter out TCGA header comment lines starting with '#'
                lines = (line for line in f if not line.startswith('#'))
                reader = csv.DictReader(lines, delimiter='\t')
                
                file_count = 0
                for row in reader:
                    v_class = str(row.get('Variant_Classification', row.get('consequence', row.get('effect', '')))).lower()
                    
                    is_noncoding = any(keyword in v_class for keyword in NONCODING_KEYWORDS)

                    if is_noncoding:
                        out_row = {
                            'Chromosome': row.get('Chromosome', row.get('chr', row.get('chrom', ''))),
                            'Start_Position': row.get('Start_Position', row.get('start', row.get('pos', ''))),
                            'End_Position': row.get('End_Position', row.get('end', '')),
                            'Reference_Allele': row.get('Reference_Allele', row.get('ref', row.get('reference', ''))),
                            'Tumor_Seq_Allele2': row.get('Tumor_Seq_Allele2', row.get('alt', row.get('tumor_seq_allele1', ''))),
                            'Variant_Classification': row.get('Variant_Classification', row.get('consequence', '')),
                            'Hugo_Symbol': row.get('Hugo_Symbol', row.get('gene', row.get('gene_symbol', ''))),
                            'Tumor_Sample_Barcode': row.get('Tumor_Sample_Barcode', row.get('sample', '')),
                            'Cancer_Type': cancer_type
                        }
                        writer.writerow(out_row)
                        file_count += 1
                        total_count += 1

                print(f"    -> Extracted {file_count} non-coding variants.")

        except Exception as e:
            print(f"    [!] Error reading file ({file}): {e}")

print(f"\n[+] SUCCESS: Merged non-coding dataset saved to: {OUTPUT_FILE}")
print(f"[+] Total non-coding variants collected: {total_count}")