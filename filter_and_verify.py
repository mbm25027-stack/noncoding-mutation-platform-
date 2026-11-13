import urllib.request
import gzip
import os

url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
gz_file = "clinvar.vcf.gz"
output_file = "blind_test_set_nc.vcf"

print(f"Current working folder: {os.getcwd()}")

# Step 1: Download the file automatically if it is missing
if not os.path.exists(gz_file):
    print("\nDownloading ClinVar dataset... (this might take a few minutes)")
    urllib.request.urlretrieve(url, gz_file)
    print("Download complete!")
else:
    print("\nClinVar dataset already found on your system.")

# Step 2: Filter for pathogenic non-coding variants
keywords = [
    "regulatory_region_variant", 
    "transcript_ablation", 
    "upstream_gene_variant", 
    "non_coding_transcript_variant"
]

matched_count = 0
print(f"Filtering VCF to extract 85 non-coding variants...")

with gzip.open(gz_file, 'rt', encoding='utf-8') as f_in, open(output_file, 'w', encoding='utf-8') as f_out:
    for line in f_in:
        # Keep header lines
        if line.startswith("#"):
            f_out.write(line)
        # Search for pathogenic regulatory variants
        elif "Pathogenic" in line:
            line_lower = line.lower()
            if any(keyword in line_lower for keyword in keywords):
                f_out.write(line)
                matched_count += 1
                if matched_count >= 85:
                    break

# Step 3: Verify the output
print("\n--- VERIFICATION REPORT ---")
if os.path.exists(output_file):
    print(f"SUCCESS: The file '{output_file}' has been successfully created.")
    print(f"Total variants extracted: {matched_count}")
    
    # Check file size to ensure it isn't empty
    file_size_kb = os.path.getsize(output_file) / 1024
    print(f"File size: {file_size_kb:.2f} KB")
else:
    print("FAILED: The file was not created. Check your folder permissions.")