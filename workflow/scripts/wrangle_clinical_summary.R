# Purpose:
# 1) Clean columns so that they are easy to access for rest of analysis
# 2) Format survival data for ease of downstream analysis
#    i) OS survival stored in vital_status, os columns
#    ii) EFS survival stored in efs_status, efs columns
# 3) Subcategorized patients from 'Other' group into finer cytogenetic groups
#     t(6;9), t(3;5)(q25;q34), del5q, del7q, del9q, monosomy5/7, trisomy8/21, -Y, -X
# New groupings are contained in the column 'cytogenetic_subgroup'

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

clinical_summary <- read_csv(snakemake@input[["clinical_summary"]])

# Cleaning column names for ease of access
formatted_clinical_summary <- clinical_summary %>%
  rename_with(~str_replace_all(., c(" " = "_", "/" = "_", "\\?" = ""))) %>%
  rename_with(~tolower(.)) %>%
  rename("patient_id" = "target_usi") %>%
  rename(sex = gender) %>%
  rename(wbc = wbc_at_diagnosis)

# Formatting survival status for ease of downstream analysis
# There are two forms of survival data that can be used: Overall survival or Event free survival
# vital_status dictates alive/dead status, needs to be reformatted
# first_event dictates if patient has either progressed (Relapse, Induction failure, Death without remission, Death)
# or not (censored)
# For EFS will use a binary classification, however perhaps may be pertinent to do competing risks somewhere down the road.
formatted_clinical_summary <- formatted_clinical_summary %>%
  mutate(os_status = ifelse(str_detect(vital_status, "Alive"), 0, 1)) %>%
  mutate(efs_status = ifelse(str_detect(first_event, "Censored"), 0, 1)) %>%
  rename(
    os_time = overall_survival_time_in_days,
    efs_time = event_free_survival_time_in_days
  )

# To determine patient subgroup, have a column called primary_cytogenetic_code
# Currently, this stratifies patients into 6 subgroups (inv16, t(8;21), Normal, MLL and 'Other' being a large contingent)
# Goal is to further delineate Other group based on observed cytogenetic abnormalities
# There is another column 'cytogenetic_complexity' that denotes number of chromosomal rearrangements.
# To assign patients to new groups, if part of 'Other', will subcategorize by cytogenetic complexity
# 1 => Patient assigned to one of t(6;9), t(3;5)(q25;q34), del5q, del7q, del9q, monosomy5/7, trisomy8/21, -Y, -X
# >=3 => Complex karyotype
# 2 => Other/mixed

subtypes <- c(
  "t(6;9)",
  "t(3;5)(q25;q34)",
  "del5q",
  "del7q",
  "del9q",
  "monosomy_5",
  "monosomy_7",
  "trisomy_8",
  "trisomy_21",
  "minus_y",
  "minus_x"
)

other_patients <- formatted_clinical_summary %>%
  filter(primary_cytogenetic_code == "Other")

complex_karyotype_patients <- other_patients %>%
  filter(cytogenetic_complexity >= 3)

other_mutations <- other_patients %>%
  filter(cytogenetic_complexity < 3) %>%
  select(patient_id, all_of(subtypes)) %>%
  pivot_longer(
    !patient_id,
    names_to = "cytogenetic_subgroup",
    values_to = "status"
  ) %>%
  filter(status == "Yes")

# Patients with count greater than 1 will be Other/combined
other_combined_patients <- other_mutations %>%
  group_by(patient_id) %>%
  count() %>%
  filter(n > 1)

# These are patients that have only one of the cytogenetic abnormalities outlined in the variable subtypes
# These will now constitute their own cytogenetic catgeory
mono_patients_groups <- other_mutations %>%
  filter(!patient_id %in% other_combined_patients[["patient_id"]])

# Storing mono-grouped patients into dictionary so groupings can be accessed in case_when call
mono_patients <- mono_patients_groups[["cytogenetic_subgroup"]]
names(mono_patients) <- mono_patients_groups[["patient_id"]]

formatted_clinical_summary <- formatted_clinical_summary %>%
  mutate(
    cytogenetic_subgroup = case_when(
      patient_id %in% complex_karyotype_patients[["patient_id"]] ~ "complex_karyotype",
      patient_id %in% other_combined_patients[["patient_id"]] ~ "other/combined_abnormalities",
      patient_id %in% names(mono_patients) ~ mono_patients[patient_id],
      .default = primary_cytogenetic_code
    )
  )

write_csv(formatted_clinical_summary, snakemake@output[["formatted_clinical_summary"]])
