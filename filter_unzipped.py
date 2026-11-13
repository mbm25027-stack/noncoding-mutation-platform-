import os

input_file = "clinvar.vcf" 
output_file = "blind_test_set_nc.vcf"

# Expanded to include common non-coding regions to guarantee 85 variants
keywords = [
    "regulatory_region_variant", 
    "upstream_gene_variant", 
    "non_coding_transcript_variant",
    "intron_variant",
    "intergenic_variant",
    "3_prime_utr_variant",
    "5_prime_utr_variant",
    "splice_region_variant"
]

matched_count = 0
print("Scanning with expanded non-coding keywords... hang tight.")

try:
    with open(input_file, 'rt', encoding='utf-8') as f_in, open(output_file, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            if line.startswith("#"):
                f_out.write(line)
                continue
            
            line_lower = line.lower()
            if "pathogenic" in line_lower:
                if any(keyword in line_lower for keyword in keywords):
                    f_out.write(line)
                    matched_count += 1
                    if matched_count >= 85:
                        break

    print(f"\nSUCCESS: Created '{output_file}'")
    print(f"Total variants extracted: {matched_count}")
except FileNotFoundError:
    print("Error: Could not find 'clinvar.vcf'. Make sure it's in the same folder.")