library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggplot2)

load("LabNotebook/Chap3/ms_phyloseq.Rdata")

sample_data(ms_phyloseq)$group5 <- case_when(
  sample_data(ms_phyloseq)$treatments == "Dimethyl fumarate" ~ "DMF",
  
  sample_data(ms_phyloseq)$treatments %in% c(
    "Glatiramer acetate",
    "Interferon"
  ) ~ "GA_IFN",
  
  sample_data(ms_phyloseq)$treatments %in% c(
    "ocrevus(rituxan)",
    "Fingolimod"
  ) ~ "Ocrevus_Fingolimod",
  
  sample_data(ms_phyloseq)$treatments == "Untreated" ~ "Untreated_PMS",
  sample_data(ms_phyloseq)$treatments == "Control" ~ "Control",
  
  TRUE ~ NA_character_
)

table(sample_data(ms_phyloseq)$group5, useNA = "ifany")

phyloseq_groups <- subset_samples(ms_phyloseq, !is.na(group5))
phyloseq_groups <- prune_taxa(taxa_sums(phyloseq_groups) > 0, phyloseq_groups)
phyloseq_groups_rel <- transform_sample_counts(phyloseq_groups, function(x) x / sum(x))

dmf <- subset_samples(phyloseq_groups_rel, group5 == "DMF")
ga_ifn <- subset_samples(phyloseq_groups_rel, group5 == "GA_IFN")
ocrev_fing <- subset_samples(phyloseq_groups_rel, group5 == "Ocrevus_Fingolimod")
pms_untreated <- subset_samples(phyloseq_groups_rel, group5 == "Untreated_PMS")
control <- subset_samples(phyloseq_groups_rel, group5 == "Control")

dmf <- prune_taxa(taxa_sums(dmf) > 0, dmf)
ga_ifn <- prune_taxa(taxa_sums(ga_ifn) > 0, ga_ifn)
ocrev_fing <- prune_taxa(taxa_sums(ocrev_fing) > 0, ocrev_fing)
pms_untreated <- prune_taxa(taxa_sums(pms_untreated) > 0, pms_untreated)
control <- prune_taxa(taxa_sums(control) > 0, control)

dmf_ASVs <- core_members(dmf, detection = 0.001, prevalence = 0.5)
ga_ifn_ASVs <- core_members(ga_ifn, detection = 0.001, prevalence = 0.5)
ocrev_fing_ASVs <- core_members(ocrev_fing, detection = 0.001, prevalence = 0.5)
pms_untreated_ASVs <- core_members(pms_untreated, detection = 0.001, prevalence = 0.5)
control_ASVs <- core_members(control, detection = 0.001, prevalence = 0.5)

venn_list_2 <- (x=list(
  DMF = dmf_ASVs,
  Immunomodulator = ga_ifn_ASVs,
  Control = control_ASVs,
  T_and_B_Cell = ocrev_fing_ASVs,
  Untreated_PMS = pms_untreated_ASVs
  
))

venn_3_groups <- ggVennDiagram(venn_list_2,
                                  label = "count",
                                  label_alpha = 0
) +
  scale_fill_gradient(low = "grey90", high = "steelblue")

ggsave("LabNotebook/Chap6/venn_3_groups.png", venn_3_groups, width = 12, height =15)
