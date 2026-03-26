# load packages
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggrepel)

# load in data (use non-rarified data)
load("LabNotebook/Chap3/ms_phyloseq.RData")

ms_filtered <- subset_samples(ms_phyloseq, disease_course_control != "Control_RRMS")

#### COMPARISON 1: Treated vs Untreated (RAW COUNTS) ####

cat("\n===== DESEQ2 WITH RAW COUNTS: TREATED VS UNTREATED =====\n")

# Subset to PMS only (no controls)
ms_pms_raw <- subset_samples(ms_filtered, 
                             treatment_status %in% c("Treated", "Untreated"))

cat("Sample sizes:\n")
print(table(sample_data(ms_pms_raw)$treatment_status))

# Add pseudocount
ms_pms_plus1 <- transform_sample_counts(ms_pms_raw, function(x) x + 1)

# DESeq2
ms_deseq_pms <- phyloseq_to_deseq2(ms_pms_plus1, ~ treatment_status)
ms_deseq_pms <- estimateSizeFactors(ms_deseq_pms)
ms_deseq_pms <- estimateDispersionsGeneEst(ms_deseq_pms)
dispersions(ms_deseq_pms) <- mcols(ms_deseq_pms)$dispGeneEst
ms_results_pms <- nbinomWaldTest(ms_deseq_pms)

res_treated_untreated <- results(ms_results_pms, 
                                 tidy = TRUE,
                                 contrast = c("treatment_status", "Treated", "Untreated"))

# Results
cat("\n=== Treated vs Untreated Results (RAW COUNTS) ===\n")
cat("Total ASVs:", nrow(res_treated_untreated), "\n")
cat("Smallest padj:", min(res_treated_untreated$padj, na.rm = TRUE), "\n")
cat("ASVs with padj < 0.05:", sum(res_treated_untreated$padj < 0.05, na.rm = TRUE), "\n")
cat("ASVs with padj < 0.1:", sum(res_treated_untreated$padj < 0.1, na.rm = TRUE), "\n")

cat("\nTop 10 ASVs:\n")
res_treated_untreated %>% 
  arrange(padj) %>% 
  select(row, baseMean, log2FoldChange, pvalue, padj) %>%
  head(10) %>%
  print()


#### COMPARISON 2: Treated vs Control (RAW COUNTS) ####

cat("\n===== DESEQ2 WITH RAW COUNTS: TREATED VS CONTROL =====\n")

ms_treated_ctrl_raw <- subset_samples(ms_filtered, 
                                      treatment_status %in% c("Treated", "Control"))

cat("Sample sizes:\n")
print(table(sample_data(ms_treated_ctrl_raw)$treatment_status))

ms_treated_ctrl_plus1 <- transform_sample_counts(ms_treated_ctrl_raw, function(x) x + 1)

ms_deseq_treated_ctrl <- phyloseq_to_deseq2(ms_treated_ctrl_plus1, ~ treatment_status)
ms_deseq_treated_ctrl <- estimateSizeFactors(ms_deseq_treated_ctrl)
ms_deseq_treated_ctrl <- estimateDispersionsGeneEst(ms_deseq_treated_ctrl)
dispersions(ms_deseq_treated_ctrl) <- mcols(ms_deseq_treated_ctrl)$dispGeneEst
ms_results_treated_ctrl <- nbinomWaldTest(ms_deseq_treated_ctrl)

res_treated_ctrl_raw <- results(ms_results_treated_ctrl, 
                                tidy = TRUE,
                                contrast = c("treatment_status", "Treated", "Control"))

cat("\n=== Treated vs Control Results (RAW COUNTS) ===\n")
cat("Total ASVs:", nrow(res_treated_ctrl_raw), "\n")
cat("Smallest padj:", min(res_treated_ctrl_raw$padj, na.rm = TRUE), "\n")
cat("ASVs with padj < 0.05:", sum(res_treated_ctrl_raw$padj < 0.05, na.rm = TRUE), "\n")
cat("ASVs with padj < 0.1:", sum(res_treated_ctrl_raw$padj < 0.1, na.rm = TRUE), "\n")

cat("\nTop 10 ASVs:\n")
res_treated_ctrl_raw %>% 
  arrange(padj) %>% 
  select(row, baseMean, log2FoldChange, pvalue, padj) %>%
  head(10) %>%
  print()


#### COMPARISON 3: Untreated vs Control (RAW COUNTS) ####

cat("\n===== DESEQ2 WITH RAW COUNTS: UNTREATED VS CONTROL =====\n")

ms_untreated_ctrl_raw <- subset_samples(ms_filtered, 
                                        treatment_status %in% c("Untreated", "Control"))

cat("Sample sizes:\n")
print(table(sample_data(ms_untreated_ctrl_raw)$treatment_status))

ms_untreated_ctrl_plus1 <- transform_sample_counts(ms_untreated_ctrl_raw, function(x) x + 1)

ms_deseq_untreated_ctrl <- phyloseq_to_deseq2(ms_untreated_ctrl_plus1, ~ treatment_status)
ms_deseq_untreated_ctrl <- estimateSizeFactors(ms_deseq_untreated_ctrl)
ms_deseq_untreated_ctrl <- estimateDispersionsGeneEst(ms_deseq_untreated_ctrl)
dispersions(ms_deseq_untreated_ctrl) <- mcols(ms_deseq_untreated_ctrl)$dispGeneEst
ms_results_untreated_ctrl <- nbinomWaldTest(ms_deseq_untreated_ctrl)

res_untreated_ctrl_raw <- results(ms_results_untreated_ctrl, 
                                  tidy = TRUE,
                                  contrast = c("treatment_status", "Untreated", "Control"))

cat("\n=== Untreated vs Control Results (RAW COUNTS) ===\n")
cat("Total ASVs:", nrow(res_untreated_ctrl_raw), "\n")
cat("Smallest padj:", min(res_untreated_ctrl_raw$padj, na.rm = TRUE), "\n")
cat("ASVs with padj < 0.05:", sum(res_untreated_ctrl_raw$padj < 0.05, na.rm = TRUE), "\n")
cat("ASVs with padj < 0.1:", sum(res_untreated_ctrl_raw$padj < 0.1, na.rm = TRUE), "\n")

cat("\nTop 10 ASVs:\n")
res_untreated_ctrl_raw %>% 
  arrange(padj) %>% 
  select(row, baseMean, log2FoldChange, pvalue, padj) %>%
  head(10) %>%
  print()


#### FINAL SUMMARY ####

cat("\n========================================\n")
cat("   FINAL SUMMARY (RAW COUNTS)\n")
cat("========================================\n\n")

cat("COMPARISON 1: Treated vs Untreated\n")
cat("  Significant ASVs (padj < 0.05):", sum(res_treated_untreated$padj < 0.05, na.rm = TRUE), "\n\n")

cat("COMPARISON 2: Treated vs Control\n")
cat("  Significant ASVs (padj < 0.05):", sum(res_treated_ctrl_raw$padj < 0.05, na.rm = TRUE), "\n\n")

cat("COMPARISON 3: Untreated vs Control\n")
cat("  Significant ASVs (padj < 0.05):", sum(res_untreated_ctrl_raw$padj < 0.05, na.rm = TRUE), "\n\n")


### CREATE VOLCANO PLOTS ###
# ===== COMPARISON 1: Treated vs Untreated =====

# Filter significant ASVs
sig_treated_untreated <- res_treated_untreated %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  dplyr::rename(ASV = row)

# Volcano plot
volcano_treated_untreated <- res_treated_untreated %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 2 ~ "Large Effect (|FC| > 2)",
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Moderate Effect (|FC| > 1)",
      padj < 0.05 ~ "Small Effect",
      TRUE ~ "Not Significant"
    ),
    neg_log10_padj = -log10(padj)
  ) %>%
  filter(!is.na(padj))

plot_volcano_treated_untreated <- ggplot(volcano_treated_untreated, 
                                         aes(x = log2FoldChange, 
                                             y = neg_log10_padj, 
                                             color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(
    "Large Effect (|FC| > 2)" = "#e74c3c",
    "Moderate Effect (|FC| > 1)" = "orange",
    "Small Effect" = "#3498db",
    "Not Significant" = "gray70"
  )) +
  labs(
    title = "Differential Abundance: Treated vs Untreated PMS",
    subtitle = paste0("Significant ASVs (padj < 0.05): ", 
                      sum(res_treated_untreated$padj < 0.05, na.rm = TRUE)),
    x = "log2 Fold Change (Treated vs Untreated)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
plot_volcano_treated_untreated

ggsave("LabNotebook/Chap7/deseq2_volcano_treated_vs_untreated_RAW.png", 
       plot_volcano_treated_untreated, 
       width = 10, 
       height = 8, 
       dpi = 300)

# ===== COMPARISON 2: Treated vs Control =====

sig_treated_ctrl <- res_treated_ctrl_raw %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  dplyr::rename(ASV = row)

cat("\n=== Treated vs Control: Significant with |log2FC| > 1 ===\n")
cat("Number of ASVs:", nrow(sig_treated_ctrl), "\n")

# Volcano plot
volcano_treated_ctrl <- res_treated_ctrl_raw %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 2 ~ "Large Effect",
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Moderate Effect",
      padj < 0.05 ~ "Small Effect",
      TRUE ~ "Not Significant"
    ),
    neg_log10_padj = -log10(padj)
  ) %>%
  filter(!is.na(padj))

plot_volcano_treated_ctrl <- ggplot(volcano_treated_ctrl, 
                                    aes(x = log2FoldChange, 
                                        y = neg_log10_padj, 
                                        color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(
    "Large Effect" = "#e74c3c",
    "Moderate Effect" = "orange",
    "Small Effect" = "#3498db",
    "Not Significant" = "gray70"
  )) +
  labs(
    title = "Differential Abundance: Treated PMS vs Healthy Control",
    subtitle = paste0("Significant ASVs (padj < 0.05): ", 
                      sum(res_treated_ctrl_raw$padj < 0.05, na.rm = TRUE)),
    x = "log2 Fold Change (Treated vs Control)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
plot_volcano_treated_ctrl

ggsave("LabNotebook/Chap7/deseq2_volcano_treated_vs_control_RAW.png", 
       plot_volcano_treated_ctrl, 
       width = 10, 
       height = 8, 
       dpi = 300)

# ===== COMPARISON 3: Untreated vs Control =====

sig_untreated_ctrl <- res_untreated_ctrl_raw %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  dplyr::rename(ASV = row)

cat("\n=== Untreated vs Control: Significant with |log2FC| > 1 ===\n")
cat("Number of ASVs:", nrow(sig_untreated_ctrl), "\n")

# Volcano plot
volcano_untreated_ctrl <- res_untreated_ctrl_raw %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 2 ~ "Large Effect",
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Moderate Effect",
      padj < 0.05 ~ "Small Effect",
      TRUE ~ "Not Significant"
    ),
    neg_log10_padj = -log10(padj)
  ) %>%
  filter(!is.na(padj))

plot_volcano_untreated_ctrl <- ggplot(volcano_untreated_ctrl, 
                                      aes(x = log2FoldChange, 
                                          y = neg_log10_padj, 
                                          color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(
    "Large Effect" = "#e74c3c",
    "Moderate Effect" = "orange",
    "Small Effect" = "#3498db",
    "Not Significant" = "gray70"
  )) +
  labs(
    title = "Differential Abundance: Untreated PMS vs Healthy Control",
    subtitle = paste0("Significant ASVs (padj < 0.05): ", 
                      sum(res_untreated_ctrl_raw$padj < 0.05, na.rm = TRUE)),
    x = "log2 Fold Change (Untreated vs Control)",
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
plot_volcano_untreated_ctrl

ggsave("LabNotebook/Chap7/deseq2_volcano_untreated_vs_control_RAW.png", 
       plot_volcano_untreated_ctrl, 
       width = 10, 
       height = 8, 
       dpi = 300)