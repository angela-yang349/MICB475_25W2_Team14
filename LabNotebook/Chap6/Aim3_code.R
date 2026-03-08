library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggplot2)
library(indicspecies)

# core microbiome analysis
load("LabNotebook/Chap3/ms_phyloseq.Rdata")

# Converting to relative abundance
phyloseq_RA <- transform_sample_counts(ms_phyloseq, function(x) x / sum(x))

# Subset data into treatment and control groups
ms_control <- subset_samples(phyloseq_RA, treatments == "Control")
ms_untreated <- subset_samples(phyloseq_RA, treatments == "Untreated")
ms_treated_fingo <- subset_samples(phyloseq_RA, treatments == "Fingolimod")
ms_treated_inter <- subset_samples(phyloseq_RA, treatments == "Interferon")
ms_treated_oc <- subset_samples(phyloseq_RA, treatments == "ocrevus(rituxan)")
ms_treated_dmf <- subset_samples(phyloseq_RA, treatments == "Dimethyl fumarate")


# setting detection and prevalence thresholds
ms_control_ASVs <- core_members(ms_control, detection=0.001, prevalence = 0.5)
ms_untreated_ASVs <- core_members(ms_untreated, detection=0.001, prevalence = 0.5)
ms_treated_fingo_ASVs <- core_members(ms_treated_fingo, detection=0.001, prevalence = 0.5)
ms_treated_inter_ASVs <- core_members(ms_treated_inter, detection=0.001, prevalence = 0.5)
ms_treated_oc_ASVs <- core_members(ms_treated_oc, detection=0.001, prevalence = 0.5)
ms_treated_dmf_ASVs <- core_members(ms_treated_dmf, detection=0.001, prevalence = 0.5)

ms_venn <- ggVennDiagram(x = list(Control = ms_control_ASVs, Untreated = ms_untreated_ASVs, Fingolimod = ms_treated_fingo_ASVs, Interferon = ms_treated_inter_ASVs, Ocrevus_rituxan = ms_treated_oc_ASVs, DMF = ms_treated_dmf_ASVs ))

# save venn diagram
ggsave("LabNotebook/Chap6/ms_venn.png", ms_venn, width = 10, height = 7)
