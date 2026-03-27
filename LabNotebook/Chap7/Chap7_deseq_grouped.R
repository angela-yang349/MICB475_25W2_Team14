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


## Create Volcano Plots ##

make_volcano <- function(res_df, title_text) {
  res_df %>%
    filter(!is.na(padj)) %>%
    mutate(
      Significance = case_when(
        padj < 0.05 & abs(log2FoldChange) > 3 ~ "Large Effect (|FC| > 2.5)",
        padj < 0.05 & abs(log2FoldChange) > 1 ~ "Moderate Effect (|FC| > 1)",
        padj < 0.05 ~ "Small Effect",
        TRUE ~ "Not Significant"
      ),
      neg_log10_padj = -log10(padj),
      label = ifelse(padj < 0.05 & abs(log2FoldChange) > 3, row, NA)
    ) %>%
    ggplot(aes(x = log2FoldChange, y = neg_log10_padj, color = Significance)) +
    geom_point(alpha = 0.7) +
    geom_text_repel(aes(label = label), na.rm = TRUE, size = 3) +
    labs(
      title = title_text,
      x = "Log2 Fold Change",
      y = "-log10(adjusted p-value)"
    ) +
    theme_minimal()
}

volcano_immuno <- make_volcano(res_immuno_untreated, "Immunomodulators vs Untreated")
volcano_immuno

volcano_tb <- make_volcano(res_tb_untreated, "T/B Cell Therapies vs Untreated")
volcano_tb

ggsave("LabNotebook/Chap7/deseq2_volcano_immuno_vs_untreated.png", 
       volcano_immuno, 
       width = 10, 
       height = 8, 
       dpi = 300)

ggsave("LabNotebook/Chap7/deseq2_volcano_TBcell_vs_untreated.png", 
       volcano_tb, 
       width = 10, 
       height = 8, 
       dpi = 300)