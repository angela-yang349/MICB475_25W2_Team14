# Group 3
# Compare grouped treatments: DMF vs GA/IFN vs Ocrevus/Fingolimod

library(tidyverse)
library(phyloseq)
library(indicspecies)

# set random seed for consistent results
set.seed(100)

# load data
load("ms_rare.RData")
ms_rare_no_RRMS_ctrl <- subset_samples(ms_rare, disease_course_control != "Control_RRMS")

# glom to genus
ms_genus <- tax_glom(ms_rare_no_RRMS_ctrl, "Genus", NArm = FALSE)
ms_genus_RA <- transform_sample_counts(ms_genus, fun=function(x) x/sum(x))

# ===============================
# GROUPING: 3 groups by drug type
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

# subset samples without NA
ms_group3 <- subset_samples(ms_genus_RA, !is.na(grouped_by_drug_type))

# run ISA
isa_group3 <- multipatt(
  t(otu_table(ms_group3)), 
  cluster = sample_data(ms_group3)$grouped_by_drug_type
)

summary(isa_group3)

taxtable <- tax_table(ms_rare_no_RRMS_ctrl) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

res_group3 <- isa_group3$sign %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value < 0.05)

View(res_group3)
write.csv(res_group3, "ISA_results_group3.csv", row.names = FALSE)