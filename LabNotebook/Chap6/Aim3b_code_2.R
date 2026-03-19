# install package
install.packages("indicspecies")

# !/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

# set seed to account for randomness
set.seed(100)

# load data
load("ms_rare.RData")
ms_rare_no_RRMS_ctrl <- subset_samples(ms_rare, disease_course_control != "Control_RRMS")

# Indicator Species/Taxa Analysis
# glom to genus
ms_genus <- tax_glom(ms_rare_no_RRMS_ctrl, "Genus", NArm = FALSE)
ms_genus_RA <- transform_sample_counts(ms_genus, fun=function(x) x/sum(x))

# ===============================
# CREATE NEW GROUPINGS
# ===============================

sample_data(ms_genus_RA)$grouped_by_drug_type <- case_when(
  sample_data(ms_genus_RA)$treatments == "Dimethyl fumarate" ~ "DMF",
  
  sample_data(ms_genus_RA)$treatments %in% c(
    "Glatiramer acetate",
    "Interferon"
  ) ~ "GA_IFN",
  
  sample_data(ms_genus_RA)$treatments %in% c(
    "ocrevus(rituxan)",
    "Fingolimod"
  ) ~ "Ocrevus_Fingolimod",
  
  TRUE ~ NA_character_
)

sample_data(ms_genus_RA)$grouped_by_mechanism <- case_when(
  sample_data(ms_genus_RA)$treatments %in% c(
    "Glatiramer acetate",
    "Interferon",
    "Dimethyl fumarate"
  ) ~ "Immunomodulators",
  
  sample_data(ms_genus_RA)$treatments %in% c(
    "ocrevus(rituxan)",
    "Fingolimod"
  ) ~ "Lymphocyte",
  
  sample_data(ms_genus_RA)$treatments == "Untreated" ~ "Untreated PMS",
  sample_data(ms_genus_RA)$treatments == "Control" ~ "Healthy Control",
  
  TRUE ~ NA_character_
)

# ===============================
# ORIGINAL ISA (treatments)
# ===============================

isa_ms <- multipatt(
  t(otu_table(ms_genus_RA)), 
  cluster = sample_data(ms_genus_RA)$treatments
)

summary(isa_ms)

taxtable <- tax_table(ms_rare_no_RRMS_ctrl) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

res <- isa_ms$sign %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value < 0.05)

View(res)
write.csv(res, "ISA_results_treatments.csv", row.names = FALSE)

# ===============================
# ISA for grouped_by_drug_type
# ===============================

ms_drug <- subset_samples(ms_genus_RA, !is.na(grouped_by_drug_type))

isa_drug <- multipatt(
  t(otu_table(ms_drug)), 
  cluster = sample_data(ms_drug)$grouped_by_drug_type
)

summary(isa_drug)

res_drug <- isa_drug$sign %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value < 0.05)

View(res_drug)
write.csv(res_drug, "ISA_results_grouped_by_drug_type.csv", row.names = FALSE)

# ===============================
# ISA for grouped_by_mechanism
# ===============================

ms_mech <- subset_samples(ms_genus_RA, !is.na(grouped_by_mechanism))

isa_mech <- multipatt(
  t(otu_table(ms_mech)), 
  cluster = sample_data(ms_mech)$grouped_by_mechanism
)

summary(isa_mech)

res_mech <- isa_mech$sign %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value < 0.05)

View(res_mech)
write.csv(res_mech, "ISA_results_grouped_by_mechanism.csv", row.names = FALSE)