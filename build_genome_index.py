import os
import glob

data_dir = "data"
output_combined = os.path.join(data_dir, "encode_cCREs_GRCh38_combined.bed")

ccre_files = [
    os.path.join(data_dir, "GRCh38-cCREs.PLS.bed"),
    os.path.join(data_dir, "GRCh38-cCREs.ELS.bed"),
    os.path.join(data_dir, "GRCh38-cCREs.CA-CTCF.bed")
]

records = []

print("[+] Reading ENCODE GRCh38 regulatory BED files...")
for filepath in ccre_files:
    if os.path.exists(filepath):
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    chrom = parts[0]
                    try:
                        start = int(parts[1])
                        end = int(parts[2])
                        label = parts[3] if len(parts) > 3 else "cCRE"
                        records.append((chrom, start, end, label))
                    except ValueError:
                        continue

print(f"[+] Total elements extracted: {len(records):,}")

# Sort naturally by chromosome and start coordinate
print("[+] Sorting genomic coordinates...")
chrom_order = {f"chr{i}": i for i in range(1, 23)}
chrom_order["chrX"] = 23
chrom_order["chrY"] = 24
chrom_order["chrM"] = 25

def sort_key(item):
    chrom, start, end, _ = item
    return (chrom_order.get(chrom, 99), chrom, start, end)

records.sort(key=sort_key)

print(f"[+] Writing sorted output to {output_combined}...")
with open(output_combined, "w", encoding="utf-8") as f:
    for chrom, start, end, label in records:
        f.write(f"{chrom}\t{start}\t{end}\t{label}\n")

print("[✔] Done! Combined sorted BED file generated successfully.")