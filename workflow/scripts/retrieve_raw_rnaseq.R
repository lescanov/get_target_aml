# Purpose:
# Acquire raw target AML RNA counts and wrangle for downstream use

library(TCGAbiolinks)
library(dplyr)
library(stringr)
library(readr)

formatted_clinical_summary <- read_csv(snakemake@input[["formatted_clinical_summary"]])

# Ran into error about running out of memory
# Try this to fix it
Sys.setenv('R_MAX_VSIZE'=32000000000)

target_rnaseq <- GDCquery(
  project = "TARGET-AML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)

GDCdownload(target_rnaseq)

# There are multiple cases per patient and these must be removed before using GDCprepare
# In the case that a patient has multiple "cases", I will retain the one that first appears
# Failing to do so will cause an erro with GDCprepare
target_rnaseq$results[[1]] <- target_rnaseq$results[[1]] %>%
  filter(duplicated(cases) == FALSE)
target_rnaseq_counts <- GDCprepare(target_rnaseq, summarizedExperiment = FALSE)

# Selecting only raw reads
target_raw_rna_counts <- target_rnaseq_counts %>%
  select(gene_id, gene_name, gene_type, contains("unstranded")) %>%
  select(!contains(c("tpm", "fpkm")))

# Tidying sample barcodes and ENSEMBL IDs (removing version number)
target_raw_rna_counts <- target_raw_rna_counts %>%
  rename_with(~ str_replace_all(., c("unstranded_" = ""))) 
target_raw_rna_counts$gene_id <- gsub("\\..*", "", target_raw_rna_counts$gene_id)

# Previously ran into errors where there were multiples of a single gene in dataset
# Will retain either first occurence of gene or unique gene ids
target_raw_rna_counts <- target_raw_rna_counts %>%
  filter(duplicated(gene_id) == FALSE)

# There are sequencing metrics at the last row of the dataframe that must be pruned
# As precaution, will also filter for "gene ids" that contain the word ENSEMBL
target_raw_rna_counts <- target_raw_rna_counts %>%
  filter(!gene_id %in% c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped")) %>%
  filter(str_detect(gene_id, "ENSG")) %>%
  rename_with( ~ str_replace_all(., c("unstranded_" = "")))

write_csv(target_raw_rna_counts, snakemake@output[["target_raw_rna_counts"]])

# Patient barcode anatomy for TARGET dataset
# TARGET-20-PARFVC-09A-01R
# > 20 means AML, 21 would be AML Induction Failed
# > PARFVC - patient id
# 09A - 09 is bone marrow
# 03 is for peripheral blood
# 04 is recurrent sample bone marrow
# 10 is normal blood, 14 is normal BM
# 40 => recurrent PB
# 41 => post treatment BM
# 42 => post treatment PB
#   >A, B, C... These represent successive aliquots
# > 01R - RNA isolation from frozen tissue, O1S is from FPE embedded

# Filtering for only first replicate of RNA extracted from bone marrow
# Choosing only first replicate so that downstream statistical can be performed (assumption of independence)
# Can consider averaging counts from multiple replicates, however edgeR requires integer counts.
# Then formatting patient samples IDs to match clinical summary by trimming barcodes.
formatted_raw_rna_counts <- target_raw_rna_counts %>%
  select(gene_id, contains("09") & contains("-01R")) %>%
  rename(ensembl_id = gene_id) %>%
  rename_with(~ str_replace_all(., c("-...-..." = "")))

write_csv(formatted_raw_rna_counts, snakemake@output[["formatted_raw_rna_counts"]])
