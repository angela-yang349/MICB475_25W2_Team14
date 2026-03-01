# Load the packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggsignif)

# Load rarefied data
load("LabNotebook/Chap3/ms_rare.RData")

#### AIM 2: Alpha diversity across treatment types (with controls) ####
treatment_types_shannon_plot <-plot_richness(ms_rare, 
                                           x = "treatments", 
                                           measures = c("Shannon")) +
  xlab("Treatment Type") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  geom_point() +
  ggtitle("Alpha Diversity Across Treatment Types and Healthy Controls")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
treatment_types_shannon_plot

#save the file as a .png on your local computer
#ggsave(filename = "LabNotebook/Chap5/treatment_types_shannon_plot.png",
       treatment_types_shannon_plot,
       height=5, width=7, dpi = 300)

## Statistical test - Kruskal-Wallis rank sum test

# Extract richness metrics and metadata
Aim2_richness_estimate <- estimate_richness(ms_rare, 
                                            measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
Aim2_alpha_samp_dat <- sample_data(ms_rare)
Aim2_alpha_samp_and_richness <- data.frame(Aim2_alpha_samp_dat, Aim2_richness_estimate)

view(Aim2_alpha_samp_and_richness)

#run Kruskal-Wallis rank sum test
kruskal_treatments <- kruskal.test( Shannon ~ treatments, data = Aim2_alpha_samp_and_richness)
kruskal_treatments

#none of these differences are significant since p=0.18 (not less than 0.05), thus further analysis is not needed


#### AIM 2: Beta diversity across treatment types (with controls) ####

treatment_types_wunifrac_dist <- distance(ms_rare, method = "wunifrac")

treatment_types_wunifrac_pcoa <- ordinate(ms_rare, 
                                          method = "PCoA", 
                                          distance = treatment_types_wunifrac_dist)

# Create PCoA plot WITH healthy controls (colored by treatment type)
treatment_types_beta_plot_all <- plot_ordination(ms_rare, 
                                                 treatment_types_wunifrac_pcoa, 
                                                 color = "treatments") +
  labs(color = "Treatment Type",
       title = "Beta Diversity Across Treatment Types and Healthy Controls",
       subtitle = "PCoA with Weighted UniFrac Distance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))
treatment_types_beta_plot_all

#ggsave("LabNotebook/Chap5/treatment_types_wunifrac_pcoa_with_controls.png", 
       treatment_types_beta_plot_all, 
       width = 9, 
       height = 6, 
       dpi = 300)

# Create subset WITHOUT healthy controls
ms_rare_pms_treatments <- subset_samples(ms_rare, disease_course != "Control")

pms_treatments_wunifrac_dist <- distance(ms_rare_pms_treatments, method = "wunifrac")

pms_treatments_wunifrac_pcoa <- ordinate(ms_rare_pms_treatments, 
                                         method = "PCoA", 
                                         distance = pms_treatments_wunifrac_dist)

# Create PCoA plot WITHOUT healthy controls
treatment_types_beta_plot_pms <- plot_ordination(ms_rare_pms_treatments, 
                                                 pms_treatments_wunifrac_pcoa, 
                                                 color = "treatments") +
  labs(color = "Treatment Type",
       title = "Beta Diversity Across PMS Treatment Types (Controls Excluded)",
       subtitle = "PCoA with Weighted UniFrac Distance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 9))
treatment_types_beta_plot_pms

#ggsave("LabNotebook/Chap5/treatment_types_wunifrac_pcoa_pms_only.png", 
       treatment_types_beta_plot_pms, 
       width = 9, 
       height = 6, 
       dpi = 300)

# PERMANOVA test (treatment type effect)
treatment_types_metadata_df <- data.frame(sample_data(ms_rare))

set.seed(123)
treatment_types_permanova <- adonis2(treatment_types_wunifrac_dist ~ treatments, 
                                     data = treatment_types_metadata_df,
                                     permutations = 999)
cat("\n=== PERMANOVA: Beta Diversity Across Treatment Types ===\n")
print(treatment_types_permanova)

#differences are significant (p=0.002) so do pairwise PERMANOVA comparisons

if(treatment_types_permanova$`Pr(>F)`[1] < 0.05) {
  
  cat("\n=== Running Pairwise PERMANOVA Comparisons ===\n")
  cat("Strategy: Each treatment vs Untreated, Untreated vs Control, Each treatment vs Control\n\n")
  
  # Get list of all treatment types (excluding Untreated and Control)
  all_treatments <- unique(as.character(sample_data(ms_rare)$treatments))
  treatment_types_only <- all_treatments[!all_treatments %in% c("Untreated", "Control")]
  
  pairwise_comparisons <- data.frame(
    Comparison = character(),
    R2 = numeric(),
    F_statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  cat("=== Each Treatment Type vs Untreated ===\n")
  
  for(treatment in treatment_types_only) {
    # Skip if treatment is "Healthy Control"
    if(treatment == "Healthy Control") next
    
    # Subset to treatment vs untreated
    ms_rare_comparison <- subset_samples(ms_rare, 
                                         treatments %in% c(treatment, "Untreated"))
    
    # Check if we have samples
    if(nsamples(ms_rare_comparison) < 3) {
      cat(treatment, "vs Untreated: Insufficient samples, skipping\n")
      next
    }
    
    # Calculate distance
    comparison_dist <- distance(ms_rare_comparison, method = "wunifrac")
    comparison_metadata <- data.frame(sample_data(ms_rare_comparison))
    
    # Run PERMANOVA
    comparison_result <- adonis2(comparison_dist ~ treatments,
                                 data = comparison_metadata,
                                 permutations = 999)
    
    # Store results
    pairwise_comparisons <- rbind(pairwise_comparisons, 
                                  data.frame(
                                    Comparison = paste(treatment, "vs Untreated"),
                                    R2 = comparison_result$R2[1],
                                    F_statistic = comparison_result$F[1],
                                    p_value = comparison_result$`Pr(>F)`[1]
                                  ))
    
    cat(treatment, "vs Untreated: R² =", round(comparison_result$R2[1], 4), 
        ", p =", comparison_result$`Pr(>F)`[1], "\n")
  }
  
  cat("=== Untreated vs Healthy Control ===\n")
  
  ms_rare_untreated_ctrl <- subset_samples(ms_rare, 
                                           treatments %in% c("Untreated", "Control"))
  
  untreated_ctrl_dist <- distance(ms_rare_untreated_ctrl, method = "wunifrac")
  untreated_ctrl_metadata <- data.frame(sample_data(ms_rare_untreated_ctrl))
  
  untreated_vs_ctrl <- adonis2(untreated_ctrl_dist ~ treatments,
                               data = untreated_ctrl_metadata,
                               permutations = 999)
  
  pairwise_comparisons <- rbind(pairwise_comparisons,
                                data.frame(
                                  Comparison = "Untreated vs Healthy Control",
                                  R2 = untreated_vs_ctrl$R2[1],
                                  F_statistic = untreated_vs_ctrl$F[1],
                                  p_value = untreated_vs_ctrl$`Pr(>F)`[1]
                                ))
  
  cat("Untreated vs Healthy Control: R² =", round(untreated_vs_ctrl$R2[1], 4),
      ", p =", untreated_vs_ctrl$`Pr(>F)`[1], "\n\n")
  
  cat("=== Each Treatment Type vs Healthy Control ===\n")
  
  for(treatment in treatment_types_only) {
    
    # Skip if treatment is "Control" or "Untreated"
    if(treatment %in% c("Control", "Untreated")) next
    
    # Subset to treatment vs control
    ms_rare_comparison <- subset_samples(ms_rare, 
                                         treatments %in% c(treatment, "Control"))
    
    # Check if we have samples
    if(nsamples(ms_rare_comparison) < 3) {
      cat(treatment, "vs Healthy Control: Insufficient samples, skipping\n")
      next
    }
    # Calculate distance
    comparison_dist <- distance(ms_rare_comparison, method = "wunifrac")
    comparison_metadata <- data.frame(sample_data(ms_rare_comparison))
    
    # Run PERMANOVA
    comparison_result <- adonis2(comparison_dist ~ treatments,
                                 data = comparison_metadata,
                                 permutations = 999)
    
    # Store results
    pairwise_comparisons <- rbind(pairwise_comparisons,
                                  data.frame(
                                    Comparison = paste(treatment, "vs Healthy Control"),
                                    R2 = comparison_result$R2[1],
                                    F_statistic = comparison_result$F[1],
                                    p_value = comparison_result$`Pr(>F)`[1]
                                  ))
    
    cat(treatment, "vs Healthy Control: R² =", round(comparison_result$R2[1], 4),
        ", p =", comparison_result$`Pr(>F)`[1], "\n")
  }
  
  # FDR correction
  pairwise_comparisons$p_fdr <- p.adjust(pairwise_comparisons$p_value, 
                                         method = "BH")
  
  # Add significance indicators (FDR < 0.05)
  pairwise_comparisons$significant_fdr <- ifelse(pairwise_comparisons$p_fdr < 0.05, 
                                                 "Yes", "No")
  
  # Sort by p-value
  pairwise_comparisons <- pairwise_comparisons[order(pairwise_comparisons$p_value), ]
  
  # Display summary
  cat("\n=== Pairwise PERMANOVA Summary ===\n")
  print(pairwise_comparisons, row.names = FALSE)
  
  # Count significant comparisons
  n_sig_unadjusted <- sum(pairwise_comparisons$p_value < 0.05)
  n_sig_fdr <- sum(pairwise_comparisons$p_fdr < 0.05)
  
  # Visual summary of significant comparisons (FDR adjusted)
  if(n_sig_fdr > 0) {
    cat("\n=== Significant Comparisons (FDR < 0.05) ===\n")
    sig_comparisons <- pairwise_comparisons[pairwise_comparisons$p_fdr < 0.05, ]
    for(i in 1:nrow(sig_comparisons)) {
      cat(sig_comparisons$Comparison[i], 
          ": R² =", round(sig_comparisons$R2[i], 4),
          ", p_FDR =", round(sig_comparisons$p_fdr[i], 4), "\n")
    }
  } else {
    cat("\n=== No comparisons remained significant after FDR correction ===\n")
  }
  
} else {
  cat("\n=== Overall PMS-only PERMANOVA not significant (p ≥ 0.05) ===\n")
  cat("Pairwise comparisons not performed.\n")
}

#write.csv(pairwise_comparisons,
          "LabNotebook/Chap5/treatment_types_pairwise_permanova.csv",
          row.names = FALSE)
