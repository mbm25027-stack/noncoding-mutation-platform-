import csv
import gzip
import os
import glob

maf_file = "October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS"
output_promoters = "liver_TERT_WDR74_promoters.tsv"

liver_keywords = {"linc-jp", "lica-fr", "liver-hcc", "lirad-fr"}
liver_donors_set = set()

def open_smart(filepath):
    with open(filepath, 'rb') as f:
        is_gz = (f.read(2) == b'\x1f\x8b')
    return gzip.open(filepath, 'rt', encoding='utf-8', errors='ignore') if is_gz else open(filepath, 'rt', encoding='utf-8', errors='ignore')

# 1. Locate Metadata File
candidates = glob.glob("pcawg_specimen_histology*")
valid_files = [f for f in candidates if os.path.getsize(f) > 50000]

if not valid_files:
    print("[!] Error: No valid metadata file (> 50 KB) found.")
    exit()

metadata_file = valid_files[0]
print(f"Using metadata file: {metadata_file}")

# 2. Extract Donor/Sample IDs (Stored as Lowercase Set)
print("Extracting Liver Cancer donor and specimen IDs...")

with open_smart(metadata_file) as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    
    for row in reader:
        if not row:
            continue
        row_str_lower = " ".join(row).lower()
        if any(kw in row_str_lower for kw in liver_keywords):
            for cell in row:
                cell_clean = cell.strip().lower()
                if cell_clean and cell_clean != 'null' and len(cell_clean) > 2:
                    liver_donors_set.add(cell_clean)

print(f"Collected {len(liver_donors_set)} Liver Cancer donor/specimen ID tokens.")

# 3. Filter MAF strictly for Liver Donors
print("Filtering MAF dataset strictly for Liver Cancer promoters...")

tert_count = 0
wdr74_count = 0
processed_liver_variants = 0
sample_id_preview = []

with open_smart(maf_file) as infile, \
     open(output_promoters, 'w', encoding='utf-8', newline='') as out_prom:
    
    reader = csv.reader(infile, delimiter='\t')
    writer_prom = csv.writer(out_prom, delimiter='\t')
    
    header = next(reader)
    writer_prom.writerow(header)
    
    for i, row in enumerate(reader):
        if not row or len(row) < 8:
            continue
        
        sample_id = row[0].strip()
        sample_id_lower = sample_id.lower()
        
        if len(sample_id_preview) < 5:
            sample_id_preview.append(sample_id)
        
        # Bidirectional case-insensitive match
        is_liver = (sample_id_lower in liver_donors_set) or \
                   any(token in sample_id_lower or sample_id_lower in token for token in liver_donors_set if len(token) > 4)
        
        if not is_liver:
            continue
            
        processed_liver_variants += 1
        chrom = str(row[1]).replace('chr', '')
        
        try:
            pos = int(row[2])
            gene = row[6]
            effect = row[7]
            
            # TERT Promoter Region (chr5: 1,295,000 - 1,295,500 or 5' annotation)
            if (chrom == '5' and 1295000 <= pos <= 1295500) or (gene == 'TERT' and '5' in effect):
                writer_prom.writerow(row)
                tert_count += 1
                
            # WDR74 Promoter Region (chr11: 62,608,000 - 62,610,000 or 5' annotation)
            elif (chrom == '11' and 62608000 <= pos <= 6261000) or (gene == 'WDR74' and '5' in effect):
                writer_prom.writerow(row)
                wdr74_count += 1

        except (ValueError, IndexError):
            continue

print("\n--- Day 1 Complete (Strict Liver Cohort) ---")
print(f"Total Liver Variants Scanned: {processed_liver_variants}")
print(f"Saved dataset: {output_promoters}")
print(f"Strict Liver TERT promoter mutations: {tert_count}")
print(f"Strict Liver WDR74 promoter mutations: {wdr74_count}")

if processed_liver_variants == 0:
    print("\n[DEBUG DIAGNOSTICS]:")
    print("Sample IDs from MAF file (First 5):", sample_id_preview)
    print("Metadata tokens sample (First 5):", list(liver_donors_set)[:5])