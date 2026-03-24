# Load libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggpicrust2)
library(pheatmap)

# Set seed for reproducibility
set.seed(1)

# Load in phyloseq object
load("LabNotebook/Chap3/ms_phyloseq.Rdata") 

# remove RRMS samples
ms_phyloseq <- subset_samples(ms_phyloseq, disease_course_control != "Control_RRMS")

# Subset samples by treatment_status
ms_healthy <- subset_samples(ms_phyloseq, treatment_status == "Control")
ms_untreated <- subset_samples(ms_phyloseq, treatment_status == "Untreated")
ms_treated <- subset_samples(ms_phyloseq, treatment_status == "Treated")

# Prune taxa with zero counts
ms_healthy <- prune_taxa(taxa_sums(ms_healthy) > 0, ms_healthy)
ms_untreated <- prune_taxa(taxa_sums(ms_untreated) > 0, ms_untreated)
ms_treated <- prune_taxa(taxa_sums(ms_treated) > 0, ms_treated)

# Transform to relative abundance
ms_healthy_RA <- transform_sample_counts(ms_healthy, function(x) x / sum(x))
ms_untreated_RA <- transform_sample_counts(ms_untreated, function(x) x / sum(x))
ms_treated_RA <- transform_sample_counts(ms_treated, function(x) x / sum(x))

# Load PICRUSt2 KO/metagenome output
ko <- read.delim("data_processed/picrust/team03_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", row.names = 1)

# Subset KO matrix by sample groups
ko_healthy <- ko[, colnames(ko) %in% sample_names(ms_healthy)]
ko_untreated <- ko[, colnames(ko) %in% sample_names(ms_untreated)]
ko_treated <- ko[, colnames(ko) %in% sample_names(ms_treated)]

# Combine metadata for DAA
meta <- sample_data(ms_phyloseq) |> data.frame() |> rownames_to_column("sample_name")

# Run DAA for treatment_status (Healthy vs Untreated vs Treated)
daa_treatment <- pathway_daa(
  abundance = ko[, colnames(ko) %in% meta$sample_name],
  metadata = meta,
  group = "treatment_status",
  daa_method = "edgeR"
)

# Annotate pathways
daa_treatment_annot <- pathway_annotation(
  pathway = "KO",
  daa_results_df = daa_treatment,
  ko_to_kegg = TRUE
)

# Filter significant pathways for plotting
sig_pathways <- daa_treatment_annot |> 
  filter(p_adjust <= 0.01, abs(log2FoldChange) > 4, !is.na(pathway_name))

# Clean pathway names and classes
sig_pathways_clean <- sig_pathways |> 
  mutate(
    pathway_class = pathway_class |> str_split(";") |> 
      lapply(function(x) { paste(trimws(x[-1]), collapse = "; ") }) |> unlist(),
    pathway_name = gsub("\\s*\\[EC:[^]]+\\]$", "", pathway_name),
    pathway_class = gsub("^Unclassified;\\s*|^Not Included in Pathway or Brite;\\s*", "", pathway_class)
  )

# Barplot of log2FC
ggplot(sig_pathways_clean, aes(
  x = reorder(pathway_class, log2FoldChange), 
  y = log2FoldChange, 
  fill = log2FoldChange > 0
)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("#CC3C82", "#699CCC"), labels = c("Under-represented","Over-represented")) +
  labs(title = "Treatment Status Contrast", x = NULL, y = "log2FC") +
  theme_classic() +
  theme(legend.position = "none", axis.text.y = element_text(size = 12))

