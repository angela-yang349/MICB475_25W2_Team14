# load packages
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggrepel)

# load in data (use non-rarified data)
load("LabNotebook/Chap3/ms_phyloseq.RData")

ms_filtered <- subset_samples(ms_phyloseq, disease_course_control != "Control_RRMS")

## ADD TREATMENT MECHANISM GROUPING ##
sample_df <- data.frame(sample_data(ms_filtered))

# Create 4-group mechanism classification
sample_df$treatment_mechanism_4grp <- case_when(
  sample_df$disease_course == "Control" ~ "Healthy Control",
  sample_df$treatments == "Untreated" ~ "Untreated PMS",
  sample_df$treatments %in% c("Glatiramer acetate", "Interferon", "Dimethyl fumarate") ~ "Immunomodulators",
  sample_df$treatments %in% c("ocrevus(rituxan)", "Fingolimod") ~ "T/B Cell Therapies",
  TRUE ~ NA_character_
)

# Update phyloseq object
sample_data(ms_filtered)$treatment_mechanism_4grp <- sample_df$treatment_mechanism_4grp


## COMPARISON: Immunomodulators vs Untreated ##

ms_immuno_untreated <- subset_samples(ms_filtered,
                                      treatment_mechanism_4grp %in% c("Immunomodulators", "Untreated PMS")
)

# set reference level
sample_data(ms_immuno_untreated)$treatment_mechanism_4grp <-
  relevel(factor(sample_data(ms_immuno_untreated)$treatment_mechanism_4grp),
          ref = "Untreated PMS")

ms_immuno_untreated_plus1 <- transform_sample_counts(ms_immuno_untreated, function(x) x + 1)

dds_immuno_untreated <- phyloseq_to_deseq2(ms_immuno_untreated_plus1, ~ treatment_mechanism_4grp)
dds_immuno_untreated <- DESeq(dds_immuno_untreated)

res_immuno_untreated <- results(dds_immuno_untreated,
                                contrast = c("treatment_mechanism_4grp", "Immunomodulators", "Untreated PMS"),
                                tidy = TRUE
)

## COMPARISON: T/B cell vs Untreated ##

ms_tb_untreated <- subset_samples(ms_filtered,
                                  treatment_mechanism_4grp %in% c("T/B Cell Therapies", "Untreated PMS")
)

# set reference level
sample_data(ms_tb_untreated)$treatment_mechanism_4grp <-
  relevel(factor(sample_data(ms_tb_untreated)$treatment_mechanism_4grp),
          ref = "Untreated PMS")

ms_tb_untreated_plus1 <- transform_sample_counts(ms_tb_untreated, function(x) x + 1)

dds_tb_untreated <- phyloseq_to_deseq2(ms_tb_untreated_plus1, ~ treatment_mechanism_4grp)
dds_tb_untreated <- DESeq(dds_tb_untreated)

res_tb_untreated <- results(dds_tb_untreated,
                            contrast = c("treatment_mechanism_4grp", "T/B Cell Therapies", "Untreated PMS"),
                            tidy = TRUE
)

# function to attach taxonomy (Species preferred, fallback to Genus then Family)
add_taxonomy <- function(res_df, physeq_obj) {
  tax_df <- as.data.frame(tax_table(physeq_obj)) %>%
    rownames_to_column("row")
  
  res_df %>%
    left_join(tax_df, by = "row") %>%
    mutate(
      TaxaLabel = case_when(
        !is.na(Species) & Species != "" ~ Species,
        !is.na(Genus) & Genus != "" ~ Genus,
        !is.na(Family) & Family != "" ~ Family,
        TRUE ~ row
      )
    )
}

res_immuno_untreated_tax <- add_taxonomy(res_immuno_untreated, ms_immuno_untreated)
res_tb_untreated_tax     <- add_taxonomy(res_tb_untreated, ms_tb_untreated)

## Create Volcano Plots ##

#Volcano plot with fc cutoff of 4
make_volcano_refined <- function(res_df, title_text, top_n_labels = 10) {
  
  fc_cutoff <- 4
  padj_cutoff <- 1e-4
  
  large_label <- paste0("Large Effect (|log2FC| > ", fc_cutoff, ")")
  moderate_label <- "Moderate Effect (|log2FC| > 1)"
  small_label <- "Small Effect"
  
  df <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      Significance = case_when(
        padj < padj_cutoff & abs(log2FoldChange) > fc_cutoff ~ large_label,
        padj < 0.05 & abs(log2FoldChange) > 1 ~ moderate_label,
        padj < 0.05 ~ small_label,
        TRUE ~ "Not Significant"
      )
    )
  
  label_df <- df %>%
    filter(Significance == large_label) %>%
    arrange(padj) %>%
    slice_head(n = top_n_labels)
  
  ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text_repel(
      data = label_df,
      aes(label = TaxaLabel),
      size = 3,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = c(
      setNames("red", large_label),
      setNames("orange", moderate_label),
      setNames("blue", small_label),
      "Not Significant" = "grey"
    )) +
    labs(
      title = title_text,
      x = "Log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Significance"
    ) +
    theme_minimal()
}

volcano_immuno <- make_volcano_refined(
  res_immuno_untreated_tax,
  "Immunomodulators vs Untreated",
  top_n_labels = 20
)
volcano_immuno

volcano_tb <- make_volcano_refined(
  res_tb_untreated_tax,
  "T/B Cell Therapies vs Untreated",
  top_n_labels = 20
)
volcano_tb


#Volcano plot with fc cutoff of 3
make_volcano_strict <- function(res_df, title_text, top_n_labels = 10) {
  
  fc_cutoff <- 3
  padj_cutoff <- 1e-4
  
  large_label <- paste0("Large Effect (|log2FC| > ", fc_cutoff, ")")
  moderate_label <- paste0("Moderate Effect (1 < |log2FC| \u2264 ", fc_cutoff, ")")
  small_label <- paste0("Small Effect (|log2FC| \u2264 1)")
  ns_label <- "Not Significant"
  
  df <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      Significance = case_when(
        padj < padj_cutoff & abs(log2FoldChange) > fc_cutoff ~ large_label,
        padj < 0.05 & abs(log2FoldChange) > 1 ~ moderate_label,
        padj < 0.05 ~ small_label,
        TRUE ~ "Not Significant"
      )
    )
  
  df$Significance <- factor(df$Significance,
                            levels = c(large_label, moderate_label, small_label, ns_label)
  )
  
  y_max <- ceiling(max(df$neg_log10_padj, na.rm = TRUE))
  if (!is.finite(y_max)) y_max <- 10
  
  label_df <- df %>%
    filter(Significance == large_label) %>%
    arrange(padj) %>%
    slice_head(n = top_n_labels)
  
  ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), 
               linetype = "dashed", 
               color = "gray40",
               linewidth = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_cutoff), 
               linetype = "dashed", 
               color = "gray40",
               linewidth = 0.5) +
    geom_vline(xintercept = 0.0, 
               linetype = "solid", 
               color = "gray40",
               linewidth = 0.5) +
    geom_text_repel(
      data = label_df,
      aes(label = TaxaLabel),
      size = 4.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50"
    ) +
    scale_color_manual(values = c(
      setNames("red", large_label),
      setNames("orange", moderate_label),
      setNames("blue", small_label),
      setNames("gray70", ns_label)
    )) +
    scale_x_continuous(
      breaks = c(-6, -3, -1, 0, 1, 3, 6),
      minor_breaks = NULL,
      expand = expansion(mult = 0.1)
    ) +
    scale_y_continuous(
      breaks = seq(0, ceiling(max(df$neg_log10_padj)), by = 2),
      minor_breaks = NULL,
      limits = c(0, NA),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Significance"
    ) +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "right",
      plot.caption = element_text(hjust = 0, size = 9),
      plot.subtitle = element_text(size = 10),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA)
    )
}

volcano_immuno <- make_volcano_strict(res_immuno_untreated_tax, "Immunomodulators vs Untreated") +
  labs(title = "Immunomodulators vs Untreated")
volcano_immuno

volcano_tb <- make_volcano_strict(res_tb_untreated_tax, "T/B Cell Therapies vs Untreated") +
  labs(title = "T/B Cell Therapies vs Untreated")
volcano_tb


ggsave("LabNotebook/Chap7/deseq2_3cutoff_volcano_immuno_vs_untreated.png", 
       volcano_immuno, 
       width = 10, 
       height = 8, 
       dpi = 300)

ggsave("LabNotebook/Chap7/deseq2_3cutoff_volcano_TBcell_vs_untreated.png", 
       volcano_tb, 
       width = 10, 
       height = 8, 
       dpi = 300)