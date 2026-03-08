library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggplot2)

# core microbiome analysis
load("LabNotebook/Chap3/ms_phyloseq.Rdata")

# Converting to relative abundance
phyloseq_RA <- transform_sample_counts(ms_phyloseq, function(x) x / sum(x))

# Subset data into treatment and control groups for specific treatments venn diagram 
ms_control <- subset_samples(phyloseq_RA, treatments == "Control")
ms_untreated <- subset_samples(phyloseq_RA, treatments == "Untreated")
ms_treated_fingo <- subset_samples(phyloseq_RA, treatments == "Fingolimod")
ms_treated_inter <- subset_samples(phyloseq_RA, treatments == "Interferon")
ms_treated_oc <- subset_samples(phyloseq_RA, treatments == "ocrevus(rituxan)")
ms_treated_dmf <- subset_samples(phyloseq_RA, treatments == "Dimethyl fumarate")

# setting detection and prevalence thresholds for specific treatments venn diagram
ms_control_ASVs <- core_members(ms_control, detection=0.001, prevalence = 0.5)
ms_untreated_ASVs <- core_members(ms_untreated, detection=0.001, prevalence = 0.5)
ms_treated_fingo_ASVs <- core_members(ms_treated_fingo, detection=0.001, prevalence = 0.5)
ms_treated_inter_ASVs <- core_members(ms_treated_inter, detection=0.001, prevalence = 0.5)
ms_treated_oc_ASVs <- core_members(ms_treated_oc, detection=0.001, prevalence = 0.5)
ms_treated_dmf_ASVs <- core_members(ms_treated_dmf, detection=0.001, prevalence = 0.5)

# subset data for general venn diagram 2
ms_control_2 <- subset_samples(phyloseq_RA, treatment_status == "Control")
ms_untreated_2 <- subset_samples(phyloseq_RA, treatment_status == "Untreated")
ms_treated_2 <- subset_samples(phyloseq_RA, treatment_status == "Treated")

# setting detection and prevalence thresholds for general venn diagram 2
ms_control_2_ASVs <- core_members(ms_control_2, detection=0.001, prevalence = 0.5)
ms_untreated_2_ASVs <- core_members(ms_untreated_2, detection=0.001, prevalence = 0.5)
ms_treated_2_ASVs <- core_members(ms_treated_2, detection=0.001, prevalence = 0.5)

# venn diagrams
ms_venn_specific_treatments <- ggVennDiagram(x = list(Control = ms_control_ASVs, Untreated = ms_untreated_ASVs, Fingolimod = ms_treated_fingo_ASVs, Interferon = ms_treated_inter_ASVs, Ocrevus_rituxan = ms_treated_oc_ASVs, DMF = ms_treated_dmf_ASVs ))
ms_venn_general <-ggVennDiagram(x =list(Control = ms_control_2_ASVs, Untreated = ms_untreated_2_ASVs, Treated=ms_treated_2_ASVs)) 

# save venn diagram
ggsave("LabNotebook/Chap6/ms_venn_specific_treatments.png", ms_venn_specific_treatments, width = 10, height = 7)
ggsave("LabNotebook/Chap6/ms_venn_general.png", ms_venn_general, width = 10, height =7)
