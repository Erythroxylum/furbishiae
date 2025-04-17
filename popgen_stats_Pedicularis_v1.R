## R Codes for Pedicularis furbishiae conservation genetics ms

################################
## read libraries
library(adegenet)
library(hierfstat)
library(pegas)
library(vcfR)
library(dartR)
library(ggplot2)
library(reshape2)
library(readr)

################################
## Set directory
setwd("~/Downloads/")

################################
## Analyses

## read in data
vcf <- read.vcfR(file = "stacks_s69ac2mac4dp6miss40.recode.vcf")

# read imap: individual mapping to population file
pop.data.7site <- read.table("sample69_site.txt", sep = "\t", header = T)

# convert vcf to genind
genind <- vcfR2genind(vcf)

# assign individuals to sites or populations
pop(genind) <- pop.data$site

# convert to genlight
gl <- gi2gl(genind)


### Figure 2b & 2c: PCoA

#run and plot PCoA with dartR
pcoa <- gl.pcoa(gl, nfactors = 10)
#axes 1 & 2
gl.pcoa.plot(pcoa, gl, ellipse=TRUE, plevel=0.9, pop.labels='pop', axis.label.size=1, hadjust=1.5,vadjust=1, 
             xaxis = 1 , yaxis = 2)
#axes 3 & 4
gl.pcoa.plot(pcoa, gl, ellipse=TRUE, plevel=0.9, pop.labels='pop', axis.label.size=1, hadjust=1.5,vadjust=1, 
             xaxis = 3 , yaxis = 4)

# SI: axes 5 & 6
gl.pcoa.plot(pcoa, gl, ellipse=TRUE, plevel=0.9, pop.labels='pop', axis.label.size=1, hadjust=1.5,vadjust=1, 
             xaxis = 5 , yaxis = 6)
# SI: axes 7 & 8
gl.pcoa.plot(pcoa, gl, ellipse=TRUE, plevel=0.9, pop.labels='pop', axis.label.size=1, hadjust=1.5,vadjust=1, 
             xaxis = 7 , yaxis = 8)
# SI: axes 9 & 10
gl.pcoa.plot(pcoa, gl, ellipse=TRUE, plevel=0.9, pop.labels='pop', axis.label.size=1, hadjust=1.5,vadjust=1, 
             xaxis = 9 , yaxis = 10)


### Table 2: Statistics on Genetic Diversity with hierfstat

#Read vcf from stacks and <https://github.com/enormandeau/stacks_workflow/tree/master>
# read in vcf
vcf <- read.vcfR(file = "stacks_s69ac2mac4dp6miss40.recode.vcf") # full s69

# vcf to genind
genind <- vcfR2genind(vcf)

# read imap; pop assignment file
pop.data <- read.table("sample69_site_4pop_ADMIX.tsv", sep = "\t", header = T)

#assign individuals from column in map file
pop(genind) <- pop.data$pop4
pop(genind) <- pop.data$site

# convert to genlight
pop.gl <- gi2gl(genind)

#bootstrap FIS
boot <- boot.ppfis(genind, diploid = TRUE, nboot = 1000)

# write to file
output_file <- "stats_ppfis_hierfstat_site.txt"
sink(output_file)
cat("boot.ppfis call:\n")
print(boot$call)
cat("\nboot.ppfis Confidence Intervals:\n")
print(boot$fis.ci)
cat("\nLevels:\n")
print(levels(genind$pop))
sink()

#hierfstat::basic.stats
bs <- basic.stats(genind)
write.table(bs$Ho, file = "stats_Ho_hierfstat_site.txt", sep = "\t")
write.table(bs$Hs, file = "stats_Hs_hierfstat_site.txt", sep = "\t")
write.table(bs$Fis, file = "stats_Fis_hierfstat_site.txt", sep = "\t")

# Open tables and average each population column across all SNP rows.


### Table 3: bootstrap Fst with hierfstat

bootfst <- boot.ppfst(genind, diploid = TRUE, nboot = 1000)
# write to file
output_file <- "stats_ppfst_hierfstat-site.txt"
sink(output_file)
cat("bootfst call:\n")
print(bootfst$call)
cat("\nboot.ppfst out:\n")
print(bootfst)
sink()






############################################################ Other stuff not used



## pi with ambiguous removed - reduces diversity

library(ape)
library(pegas)

# Read the FASTA file
fasta_file <- "s69_sequences_for_phylogeny_first_filters_m3_p70_x0_S2.fasta"
dna <- read.dna(fasta_file, format = "fasta", as.character = TRUE)

# Count N or ambiguous IUPAC codes
ambiguous_bases <- c("r", "y", "s", "w", "k", "m", "b", "d", "h", "v", "n")
table(unlist(dna) %in% ambiguous_bases)

### pegas pi
# Calculate nucleotide diversity (π)
dna_bin <- as.DNAbin(dna)
pi_val <- nuc.div(dna_bin)

# Output the result for species:
cat("Nucleotide diversity (π):", pi_val, "\n")

# Load population map
pop_map <- read.table("sample69_site.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Match FASTA sequence names with population labels
matched_samples <- intersect(rownames(dna_bin), pop_map$sample)
dna_bin <- dna_bin[matched_samples, ]
pop_map <- pop_map[pop_map$sample %in% matched_samples, ]

# Assign population vector
pop_vector <- setNames(pop_map$site, pop_map$sample)
pop_factor <- as.factor(pop_vector[rownames(dna_bin)])

# Calculate nucleotide diversity (π) per population
pi_per_pop <- tapply(1:nrow(dna_bin), pop_factor, function(idx) {
  nuc.div(dna_bin[idx, ])
})

# Output results
print(round(pi_per_pop, 5))




############################
##### Random pseudo-phasing into haplotypes

# Load required packages
library(ape)
library(pegas)

# Step 1: Read aligned FASTA file (as character matrix)
fasta_file <- "s69_sequences_for_phylogeny_first_filters_m3_p70_x0_S2.fasta"
dna <- read.dna(fasta_file, format = "fasta", as.character = TRUE)
dna <- apply(dna, c(1, 2), toupper)  # Ensure uppercase bases
samples <- rownames(dna)
alignment_length <- ncol(dna)

# Step 2: Define IUPAC ambiguity codes (diploid)
iupac <- list(
  R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"), W = c("A", "T"),
  K = c("G", "T"), M = c("A", "C"), B = c("C", "G", "T"), D = c("A", "G", "T"),
  H = c("A", "C", "T"), V = c("A", "C", "G"), N = c("A", "C", "G", "T")
)

# Step 3: Create two pseudo-haplotypes per individual
pseudo_haps <- list()
for (i in seq_len(nrow(dna))) {
  hap1 <- hap2 <- character(alignment_length)
  for (j in seq_len(alignment_length)) {
    base <- toupper(dna[i, j])
    if (base %in% names(iupac)) {
      alleles <- sample(iupac[[toupper(base)]], 2, replace = FALSE)
      hap1[j] <- alleles[1]
      hap2[j] <- alleles[2]
    } else {
      hap1[j] <- base
      hap2[j] <- base
    }
  }
  sample_id <- samples[i]
  pseudo_haps[[paste0(sample_id, "_a")]] <- hap1
  pseudo_haps[[paste0(sample_id, "_b")]] <- hap2
}

# Step 4: Convert to DNAbin
hap_mat <- do.call(rbind, pseudo_haps)
rownames(hap_mat) <- names(pseudo_haps)
hap_bin <- as.DNAbin(hap_mat)

# Step 5: Load population map
pop_map <- read.table("sample69_site.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(pop_map) <- pop_map$sample

# Step 6: Expand population map for haplotypes
pop_map_expanded <- data.frame(
  sample = c(paste0(pop_map$sample, "_a"), paste0(pop_map$sample, "_b")),
  site = rep(pop_map$site, each = 2),
  stringsAsFactors = FALSE
)

# Step 7: Ensure consistent sample order
shared_samples <- intersect(rownames(hap_bin), pop_map_expanded$sample)
hap_bin <- hap_bin[shared_samples, ]
pop_map_expanded <- pop_map_expanded[match(shared_samples, pop_map_expanded$sample), ]

# Step 8: Assign population factors
pop_vector <- setNames(pop_map_expanded$site, pop_map_expanded$sample)
pop_factor <- factor(pop_vector[rownames(hap_bin)])

# Step 9: Define 7 populations + "Maine" composite
unique_pops <- levels(pop_factor)
maine_pops <- c("1ME-BB", "2ME-CA", "3ME-DB", "4ME-PB")
all_pops <- union(unique_pops, "Maine")

# Step 10: Calculate π per population
pi_per_pop <- sapply(all_pops, function(pop) {
  if (pop == "Maine") {
    inds <- names(pop_factor)[pop_factor %in% maine_pops]
  } else {
    inds <- names(pop_factor)[pop_factor == pop]
  }
  if (length(inds) >= 2) {
    nuc.div(hap_bin[inds, ])
  } else {
    NA
  }
})

# Step 11: Output results
cat("\nNucleotide diversity (π) per population (randomly phased):\n")
print(round(pi_per_pop, 5))

#1ME-BB  2ME-CA  3ME-DB  4ME-PB  5NB-GF  6NB-MD  8NB-BF   Maine 
#0.03004 0.02147 0.02584 0.03475 0.03207 0.02528 0.02804 0.03048 




#################################
###### bootstrap random allele phasing

library(ape)
library(pegas)

# Input files
fasta_file <- "s69_sequences_for_phylogeny_first_filters_m3_p70_x0_S2.fasta"
pop_file <- "sample69_site.txt"

# Parameters
n_iter <- 100

# Step 1: Load sequence data
dna <- read.dna(fasta_file, format = "fasta", as.character = TRUE)
dna <- apply(dna, c(1, 2), toupper)
samples <- rownames(dna)
alignment_length <- ncol(dna)

# Step 2: Load population map
pop_map <- read.table(pop_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(pop_map) <- pop_map$sample
samples <- intersect(rownames(dna), pop_map$sample)
dna <- dna[samples, ]
pop_map <- pop_map[samples, ]

# Step 3: Define IUPAC ambiguity codes
iupac <- list(
  R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"), W = c("A", "T"),
  K = c("G", "T"), M = c("A", "C"), B = c("C", "G", "T"), D = c("A", "G", "T"),
  H = c("A", "C", "T"), V = c("A", "C", "G"), N = c("A", "C", "G", "T")
)

# Step 4: Define population groups
unique_pops <- unique(pop_map$site)
maine_pops <- c("1ME-BB", "2ME-CA", "3ME-DB", "4ME-PB")
all_pops <- c(unique_pops, "Maine", "All")

# Step 5: Prepare result matrix
pi_matrix <- matrix(NA, nrow = n_iter, ncol = length(all_pops))
colnames(pi_matrix) <- all_pops

# Step 6: Iterative pseudo-phasing + π calculation
set.seed(NULL)  # Allow natural RNG

for (iter in 1:n_iter) {
  pseudo_haps <- list()
  for (i in seq_len(nrow(dna))) {
    hap1 <- hap2 <- character(alignment_length)
    for (j in seq_len(alignment_length)) {
      base <- toupper(dna[i, j])
      if (base %in% names(iupac)) {
        alleles <- sample(iupac[[base]], 2, replace = FALSE)
        hap1[j] <- alleles[1]
        hap2[j] <- alleles[2]
      } else {
        hap1[j] <- base
        hap2[j] <- base
      }
    }
    sample_id <- rownames(dna)[i]
    pseudo_haps[[paste0(sample_id, "_a")]] <- hap1
    pseudo_haps[[paste0(sample_id, "_b")]] <- hap2
  }
  
  hap_mat <- do.call(rbind, pseudo_haps)
  rownames(hap_mat) <- names(pseudo_haps)
  hap_bin <- as.DNAbin(hap_mat)
  
  # Assign populations
  pop_map_expanded <- data.frame(
    sample = c(paste0(pop_map$sample, "_a"), paste0(pop_map$sample, "_b")),
    site = rep(pop_map$site, each = 2),
    stringsAsFactors = FALSE
  )
  shared_samples <- intersect(rownames(hap_bin), pop_map_expanded$sample)
  hap_bin <- hap_bin[shared_samples, ]
  pop_map_expanded <- pop_map_expanded[match(shared_samples, pop_map_expanded$sample), ]
  pop_vector <- setNames(pop_map_expanded$site, pop_map_expanded$sample)
  pop_factor <- factor(pop_vector[rownames(hap_bin)])
  
  # π per group
  for (pop in unique_pops) {
    inds <- names(pop_factor)[pop_factor == pop]
    if (length(inds) >= 2) {
      pi_matrix[iter, pop] <- nuc.div(hap_bin[inds, ])
    }
  }
  
  # π for Maine
  maine_inds <- names(pop_factor)[pop_factor %in% maine_pops]
  if (length(maine_inds) >= 2) {
    pi_matrix[iter, "Maine"] <- nuc.div(hap_bin[maine_inds, ])
  }
  
  # π for All
  pi_matrix[iter, "All"] <- nuc.div(hap_bin)
  
  if (iter %% 50 == 0) {
    cat("Iteration", iter, "complete.\n")
  }
}

# Step 7: Summarize results
pi_summary <- apply(pi_matrix, 2, function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    lower = quantile(x, 0.025, na.rm = TRUE),
    upper = quantile(x, 0.975, na.rm = TRUE)
  )
})

pi_summary <- t(round(pi_summary, 5))
cat("\nBootstrapped π estimates (mean, SD, 95% CI):\n")
print(pi_summary)

write.csv(pi_matrix, "pi_matrix_raw.csv", row.names = FALSE)
write.csv(pi_summary, "pi_summary_bootstrapped.csv")


####### plot pi variation

desired_order <- c("All", "Maine", "1ME-BB", "2ME-CA", "3ME-DB", "4ME-PB", "5NB-GF", "6NB-MD", "8NB-BF")

pi_matrix <- pi_matrix[, desired_order]

boxplot(pi_matrix, las = 2, col = "skyblue", main = "π estimates by subpopulation", ylab = "π")


######### calculate site stats

# Read aligned DNA sequence (as character matrix)
fasta_file <- "s69_sequences_for_phylogeny_first_filters_m3_p70_x0_S2.fasta"
dna <- read.dna(fasta_file, format = "fasta", as.character = TRUE)
dna <- apply(dna, c(1, 2), toupper)

# Define standard and ambiguous bases
standard_bases <- c("A", "C", "G", "T")
ambiguous_bases <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N")

# Initialize counters
n_sites <- ncol(dna)
invariant_sites <- 0
ambiguous_sites <- 0
variable_sites <- 0
ambig_and_variable <- 0
ambig_and_invariant <- 0

# Analyze each site
for (j in 1:n_sites) {
  site_bases <- dna[, j]
  base_counts <- table(site_bases)
  
  has_ambig <- any(names(base_counts) %in% ambiguous_bases)
  standard_alleles <- intersect(names(base_counts), standard_bases)
  
  if (has_ambig) ambiguous_sites <- ambiguous_sites + 1
  
  if (length(standard_alleles) == 1 && !has_ambig) {
    invariant_sites <- invariant_sites + 1
  } else if (length(standard_alleles) > 1) {
    variable_sites <- variable_sites + 1
    if (has_ambig) ambig_and_variable <- ambig_and_variable + 1
  } else if (length(standard_alleles) == 1 && has_ambig) {
    ambig_and_invariant <- ambig_and_invariant + 1
  }
}

# Summary
cat("Total sites:", n_sites, "\n")
cat("Invariant (standard only):", invariant_sites, "\n")
cat("Ambiguous sites (any ambiguity):", ambiguous_sites, "\n")
cat("Variable (multiple standard alleles):", variable_sites, "\n")
cat("Ambiguous + Variable sites:", ambig_and_variable, "\n")
cat("Ambiguous + Invariant sites:", ambig_and_invariant, "\n")

# Define ambiguous IUPAC codes
ambiguous_bases <- c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N")

# Flatten the DNA matrix and count ambiguous characters
total_ambiguous <- sum(dna %in% ambiguous_bases)

# Output result
total_bases <- length(dna)
percent_ambiguous <- 100 * total_ambiguous / total_bases

cat("Total ambiguous nucleotides in alignment:", total_ambiguous, ", (", round(percent_ambiguous, 2), "% out of", length(dna), "total bases)")
                            