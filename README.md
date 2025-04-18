# furbishiae

Bash and R scripts for conservation genetics analysis of *Pedicularis furbishiae*.

---

### 📂 Data

- **STACKS output SNP data**  
  `first_filters_m3_p70_x0_S2.singleton.unlinked_0.5_100kbp.vcf.gz`  
  [Download](https://drive.google.com/file/d/1Y8Ew879WDDAhEZCuA0PPKnGqKwD_LInA/view?usp=sharing)

---

### 🧪 Analysis Scripts

#### Step 1: Filter VCF output from STACKS using VCFtools  
(After following our methods using [stacks_workflow](https://github.com/enormandeau/stacks_workflow))

```bash
vcftools --gzvcf first_filters_m3_p70_x0_S2.singleton.unlinked_0.5_100kbp.vcf.gz \
  --keep sample69.txt \
  --min-alleles 2 --max-alleles 2 \
  --mac 4 \
  --minDP 6 \
  --max-missing 0.5 \
  --recode
  --out stacks_s69ac2mac4dp6miss40
```

- **Or download cleaned SNP data**  
  `stacks_s69ac2mac4dp6miss40.recode.vcf.gz`  
  [Download](https://drive.google.com/file/d/15tJFhJ9qe3nL2_6aTbCCPKj_mKdQJ1Xc/view?usp=sharing)


#### Step 2: LD pruning with PLINK for ADMIXTURE

Generate .bed for ADMIXTURE:

```bash
./run_generate_admixture_files.sh file.recode.vcf
```

- **Or download LD-pruned SNP data (for ADMIXTURE)**  
  `stacks_s69ac2mac4dp6miss40_ld50kb.recode.vcf.gz`  
  [Download](https://drive.google.com/file/d/16iwUK8R8GXZu82t2a5cKkH6ncZ-uNZgT/view?usp=sharing)
  
---

#### Step 3: Run ADMIXTURE

```bash
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
```

> Plotting uses the R script `admixture_multiple_k_labels_DW.R`  
> (adapted from Joana Meier and Eric Normandeau)

---

#### Step 4: PCA and Population Genetic Stats in R

- **Download sequence data for π calculation**  
  `s69_sequences_for_phylogeny_first_filters_m3_p70_x0_S2.fasta`  
  [Download](https://drive.google.com/file/d/1GPehOBpxqyjjUoU2sPLVo7qG42fX_UUc/view?usp=sharing)

  Open "popgen_stats_Pedicularis_v1.R" and run codes.

---
