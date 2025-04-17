# furbishiae

Bash and R scripts for conservation genetics analysis of *Pedicularis furbishiae*.

---

### 📂 Data

- **STACKS output SNP data**  
  `first_filters_m3_p70_x0_S2.singleton.unlinked_0.5_100kbp.vcf.gz`  
  [Download](https://drive.google.com/file/d/1Y8Ew879WDDAhEZCuA0PPKnGqKwD_LInA/view?usp=sharing)

- **Cleaned SNP data**  
  `stacks_s69ac2mac4dp6miss40.recode.vcf.gz`  
  [Download](https://drive.google.com/file/d/15tJFhJ9qe3nL2_6aTbCCPKj_mKdQJ1Xc/view?usp=sharing)

- **LD-pruned SNP data (for ADMIXTURE)**  
  `stacks_s69ac2mac4dp6miss40_ld50kb.recode.vcf.gz`  
  [Download](https://drive.google.com/file/d/16iwUK8R8GXZu82t2a5cKkH6ncZ-uNZgT/view?usp=sharing)

- **Sequence data for π calculation**  
  `s69_sequences_for_phylogeny_first_filters_m3_p70_x0_S2.fasta`  
  [Download](https://drive.google.com/file/d/1GPehOBpxqyjjUoU2sPLVo7qG42fX_UUc/view?usp=sharing)

---

### 🧪 Analysis Scripts

#### Step 1: Filter VCF output from STACKS using VCFtools  
(based on [stacks_workflow](https://github.com/enormandeau/stacks_workflow))

```bash
vcftools --gzvcf first_filters_m3_p70_x0_S2.singleton.unlinked_0.5_100kbp.vcf.gz \
  --keep sample69.txt \
  --min-alleles 2 --max-alleles 2 \
  --mac 4 \
  --minDP 6 \
  --max-missing 0.5 \
  --out stacks_s69ac2mac4dp6miss40
```

#### Step 2: LD pruning with PLINK for ADMIXTURE

```bash
plink --vcf stacks_s69ac2mac4dp6miss40.recode.vcf \
  --indep-pairwise 50 5 0.2 \
  --out pruned
```

Run:

```bash
./run_generate_admixture_files.sh
```

> Plotting uses the R script `admixture_multiple_k_labels_DW.R`  
> (adapted from Joana Meier and Eric Normandeau)

---

#### Step 3: Run ADMIXTURE

```bash
./run_admixture.sh
```

---

#### Step 4: PCA and Population Genetic Stats in R

```r
source("popgen_stats_Pedicularis_v1.R")
```

---
