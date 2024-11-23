# Purpose
# Retrieve raw miRNA-seq counts

library(TCGAbiolinks)
library(stringr)
library(readr)
library(dplyr)

formatted_clinical_summary <- read_csv(snakemake@input[["formatted_clinical_summary"]])

target_mirna <- GDCquery(
  project = "TARGET-AML",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling"
)

GDCdownload(target_mirna, files.per.chunk = 50)
target_raw_mirna_counts <- GDCprepare(target_mirna)

target_raw_mirna_counts <- target_raw_mirna_counts %>%
  select(miRNA_ID, contains("read_count")) %>%
  rename_with(~ str_replace_all(., c("read_count_" = ""))) %>%
  rename(mirna_id = miRNA_ID)

# Removing any possible QC statistics
target_raw_mirna_counts <- target_raw_mirna_counts %>%
  filter(!mirna_id %in% c("N_ambiguous", "N_multimapping", "N_noFeature", "N_unmapped"))

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

formatted_raw_mirna_counts <- target_raw_mirna_counts %>%
  select(
    mirna_id,
    contains("TARGET-20") &
      contains("-09A") &
      contains("-01R") &
      contains(formatted_clinical_summary[["patient_id"]])
  ) %>%
  rename_with(~ str_replace_all(., c("-...-..." = "")))

write_csv(formatted_raw_mirna_counts, snakemake@output[["formatted_raw_mirna_counts"]])

