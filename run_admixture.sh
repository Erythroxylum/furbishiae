# Bash scripts for Pedicularis furbishiae population genetics analyses:

## step 1: filtering vcf output from STACKS and <https://github.com/enormandeau/stacks_workflow/tree/master>
# set mac and dp and missing 40%, remove admix
# on command line:
vcftools --vcf first_filters_m3_p70_x0_S2.singleton.unlinked_0.5_100kbp.vcf --keep sample69.txt --min-alleles  2 --max-alleles 2 --mac 4 --minDP 6 --max-missing 0.5 --out stacks_s69ac2mac4dp6miss40


# FIGURE 2a:
# Filter SNPs for LD using PLINK and produce ADMIXTURE files:

# after running run_generate_admixture_files.sh 

# run ADMIXTURE using GNU Parallel
seq 10 | parallel admixture stacks_s69ac2mac4dp6miss40_ld50kb.bed {} -j3 --cv -C 0.00001 -c 0.000000000001 \> {}.log &

## get CVerror
grep -h CV *.log | sort -V | awk '{print $3,$4}' | cut -d "=" -f 2  | perl -pe 's/\)://' | awk '{print $2,$1}' > cverror.txt

## plot ADMIXTURE barplots
# 7site with labels
rscript /Users/dawsonwhite/Library/Mobile\ Documents/com\~apple\~CloudDocs/R/admixture_multiple_k_labels_DW.R \
-p stacks_s69ac2mac4dp6miss40_ld50kb \
-i ../sample69_site.txt \
-k 10 -m 2 -l 1ME-BB,2ME-CA,3ME-DB,4ME-PB,5NB-GF,6NB-MD,8NB-BF