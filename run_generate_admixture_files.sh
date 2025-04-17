#!/bin/bash

# Function to display an error message and exit
error_exit() {
  echo "$1" 1>&2
  exit 1
}

# Check if filename is provided as a command-line argument
if [ -z "$1" ]; then
  error_exit "Error: No filename provided. Usage: $0 <filename.recode.vcf>"
fi

# Extract the base filename without .recode.vcf extension
filename=$(basename "$1" .recode.vcf)

# Check if the original file with .recode.vcf exists
if [ ! -f "${filename}.recode.vcf" ]; then
  error_exit "Error: VCF file '${filename}.recode.vcf' not found."
fi

# Create directory based on the input filename
admixture_dir="admixture_${filename}"
mkdir -p "$admixture_dir"

# Generate PLINK files from VCF
echo "Generating PLINK files from VCF..."
plink --allow-extra-chr \
      --vcf "${filename}.recode.vcf" \
      --recode \
      --out "$admixture_dir/${filename}" \
      || error_exit "Error: Failed to generate PLINK files from VCF."

# Perform LD pruning
echo "Performing LD pruning..."
plink --allow-extra-chr \
      --indep-pairwise 50 5 0.2 \
      --file "$admixture_dir/${filename}" \
      --out "$admixture_dir/${filename}_ld50kb" \
      || error_exit "Error: LD pruning failed."

# Check if pruning result file exists
if [ ! -f "$admixture_dir/${filename}_ld50kb.prune.in" ]; then
  error_exit "Error: LD-pruned SNP file not found."
fi

# Generate VCF with LD-pruned SNPs
echo "Generating VCF with LD-pruned SNPs..."
plink --allow-extra-chr \
      --file "$admixture_dir/${filename}" \
      --extract "$admixture_dir/${filename}_ld50kb.prune.in" \
      --recode vcf \
      --out "${filename}_ld50kb_doubleid" \
      || error_exit "Error: Failed to generate VCF with LD-pruned SNPs."
      
# If printing with double IDs, remove here with sed
#sed -e 's/ANB-1-1_ANB-1-1/ANB-1-1/' -e 's/ANB-1-2_ANB-1-2/ANB-1-2/'   "${filename}_ld50kb_doubleid.vcf" > "${filename}_ld50kb.vcf" \
#     || error_exit "Error: Failed to change doubleid"

# Generate PLINK files from LD-pruned VCF
echo "Generating PLINK files from LD-pruned VCF..."
plink --allow-extra-chr \
      --vcf "${filename}_ld50kb.vcf" \
      --make-bed \
      --out "$admixture_dir/${filename}_ld50kb" \
      || error_exit "Error: Failed to generate PLINK files from LD-pruned VCF."

# ADMIXTURE requires chromosome names to be 0 for non-human chromosomes
echo "Adjusting chromosome names for ADMIXTURE..."
awk '{$1="0"; print $0}' "$admixture_dir/${filename}_ld50kb.bim" > "$admixture_dir/${filename}_ld50kb.bim.tmp"
mv "$admixture_dir/${filename}_ld50kb.bim.tmp" "$admixture_dir/${filename}_ld50kb.bim"

echo "Process completed successfully."

