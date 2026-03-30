library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggplot2)

load("LabNotebook/Chap3/ms_phyloseq.Rdata")

sample_variables(ms_phyloseq)

# mechanism of action 
sample_data(ms_phyloseq)$group4 <- case_when(
  sample_data(ms_phyloseq)$treatments %in% c(
    "Glatiramer acetate",
    "Interferon",
    "Dimethyl fumarate"
  ) ~ "Immunomodulators",
  
  sample_data(ms_phyloseq)$treatments %in% c(
    "ocrevus(rituxan)",
    "Fingolimod"
  ) ~ "Lymphocyte",
  
  sample_data(ms_phyloseq)$treatments == "Untreated" ~ "Untreated PMS",
  sample_data(ms_phyloseq)$treatments == "Control" ~ "Healthy Control",
  
  TRUE ~ NA_character_
)

table(sample_data(ms_phyloseq)$group4, useNA = "ifany")

ms_groups <- subset_samples(ms_phyloseq, !is.na(group4))
ms_groups <- prune_taxa(taxa_sums(ms_groups) > 0, ms_groups)

phyloseq_RA <- transform_sample_counts(ms_groups, function(x) x / sum(x))

ms_immunomod <- subset_samples(phyloseq_RA, group4 == "Immunomodulators")
ms_lymphocyte <- subset_samples(phyloseq_RA, group4 == "Lymphocyte")
ms_pms_untreated <- subset_samples(phyloseq_RA, group4 == "Untreated PMS")
ms_hc <- subset_samples(phyloseq_RA, group4 == "Healthy Control")

ms_immunomod <- prune_taxa(taxa_sums(ms_immunomod) > 0, ms_immunomod)
ms_lymphocyte <- prune_taxa(taxa_sums(ms_lymphocyte) > 0, ms_lymphocyte)
ms_pms_untreated <- prune_taxa(taxa_sums(ms_pms_untreated) > 0, ms_pms_untreated)
ms_hc <- prune_taxa(taxa_sums(ms_hc) > 0, ms_hc)

ms_control_ASVs <- core_members(ms_hc, detection=0.001, prevalence = 0.5)
ms_untreated_ASVs <- core_members(ms_pms_untreated, detection=0.001, prevalence = 0.5)
ms_lymphocyte_ASVs <- core_members(ms_lymphocyte, detection=0.001, prevalence = 0.5)
ms_immunomod_ASVs <- core_members(ms_immunomod, detection=0.001, prevalence = 0.5)

length(ms_immunomod_ASVs)
length(ms_lymphocyte_ASVs)
length(ms_untreated_ASVs)
length(ms_control_ASVs)

venn_list <- (x= list(
  Immunomodulator = ms_immunomod_ASVs,
  T_and_B_Cell = ms_lymphocyte_ASVs,
  Untreated_PMS = ms_untreated_ASVs,
  Control = ms_control_ASVs
))


venn_mechanism_1 <- ggVennDiagram(venn_list,
                                  label = "count",
                                  label_alpha = 0
) +
  scale_fill_gradient(low = "grey90", high = "steelblue") +
  theme(
    text = element_text(size = 30),          # overall text size
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14)
  )

ggsave("LabNotebook/Chap6/venn_mechanism_1.png", venn_mechanism_1, width = 12, height =15)

# find unique ASVs
lymph_unique <- setdiff(ms_lymphocyte_ASVs, union(ms_control_ASVs, union(ms_untreated_ASVs, ms_immunomod_ASVs)))
lymph_unique

immuno_unique <- setdiff(ms_immunomod_ASVs, union(ms_control_ASVs, union(ms_untreated_ASVs, ms_lymphocyte_ASVs)))
immuno_unique

immuno_lymph_shared <- intersect (ms_immunomod_ASVs, ms_lymphocyte_ASVs)
exclude <- union (ms_control_ASVs, ms_untreated_ASVs)
lymph_immuno_only <- setdiff(immuno_lymph_shared, exclude)
lymph_immuno_only
