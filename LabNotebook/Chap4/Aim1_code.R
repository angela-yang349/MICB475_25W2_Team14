# load packages
library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(picante)
library(ggsignif)

# load in data
load("LabNotebook/Chap3/ms_phyloseq.RData")
load("LabNotebook/Chap3/ms_rare.RData")

#### AIM 1 PRELIMINARY: Beta diversity SPMS vs PPMS (no controls) ####

# Remove control samples for PMS-only analyses
ms_rare_pms_only <- subset_samples(ms_rare, disease_course != "Control")

pms_only_wunifrac_dist <- distance(ms_rare_pms_only, method = "wunifrac")

pms_only_wunifrac_pcoa <- ordinate(ms_rare_pms_only, 
                                   method = "PCoA", 
                                   distance = pms_only_wunifrac_dist)

spms_ppms_pcoa_plot <- plot_ordination(ms_rare_pms_only, 
                                       pms_only_wunifrac_pcoa, 
                                       color = "disease_course") +
  labs(color = "Disease Course",
       title = "Beta Diversity PCoA: SPMS vs PPMS") +
  theme_classic()
spms_ppms_pcoa_plot

##only need to run the save codes once each
#ggsave("LabNotebook/Chap4/spms_ppms_wunifrac_pcoa.png",
       #spms_ppms_pcoa_plot,
       #height = 4, width = 5)

# PERMANOVA test (SPMS vs PPMS)
pms_only_metadata <- data.frame(sample_data(ms_rare_pms_only))

set.seed(123)
spms_ppms_permanova <- adonis2(pms_only_wunifrac_dist ~ disease_course, 
                               data = pms_only_metadata,
                               permutations = 999)
print(spms_ppms_permanova)


#### AIM 1: Alpha diversity treated vs untreated PMS (with controls) ####
treatment_alpha_plot <- plot_richness(ms_rare, 
                                      x = "treatment_status", 
                                      measures = c("Shannon")) +
  xlab("Treatment Status") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  geom_point() +
  ggtitle("Alpha Diversity: Treated vs Untreated PMS and Controls") +
  theme_classic()
treatment_alpha_plot

#ggsave(filename = "LabNotebook/Chap4/treatment_shannon_boxplot.png",
       #treatment_alpha_plot,
       #height = 4, width = 6)

### Statistical test - Kruskal-Wallis rank sum test
#first extract info
Aim1_richness_estimate <- estimate_richness(ms_rare, 
                                            measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
Aim1_alpha_samp_dat <- sample_data(ms_rare)
Aim1_alpha_samp_and_richness <- data.frame(Aim1_alpha_samp_dat, Aim1_richness_estimate)

view(Aim1_alpha_samp_and_richness)

#run Kruskal-Wallis rank sum test
kruskal_treatment_status <- kruskal.test( Shannon ~ treatment_status, data = Aim1_alpha_samp_and_richness)
kruskal_treatment_status

#none of these differences are significant since p=0.2533 (not less than 0.05), thus further analysis is not needed



#### AIM 1: Beta diversity treated vs untreated PMS (with controls) ####
treatment_wunifrac_dist <- distance(ms_rare, method = "wunifrac")

treatment_wunifrac_pcoa <- ordinate(ms_rare, 
                                    method = "PCoA", 
                                    distance = treatment_wunifrac_dist)

treatment_beta_plot <- plot_ordination(ms_rare, 
                                       treatment_wunifrac_pcoa, 
                                       color = "treatment_status") +
  labs(color = "Treatment Status", 
       title = "Beta Diversity: Treated vs Untreated PMS vs. Healthy Controls",
       subtitle = "PCoA with Weighted UniFrac Distance") +
  theme_classic() +
  theme(legend.position = "right")
treatment_beta_plot

#ggsave("LabNotebook/Chap4/treatment_wunifrac_pcoa.png", 
       #treatment_beta_plot, 
       #width = 8, 
       #height = 6, 
       #dpi = 300)

# PERMANOVA test (treatment effect)
metadata_all_df <- data.frame(sample_data(ms_rare))

set.seed(123)
treatment_permanova <- adonis2(treatment_wunifrac_dist ~ treatment_status, 
                               data = metadata_all_df,
                               permutations = 999)
print(treatment_permanova)

# Save overall PERMANOVA results
write.csv(as.data.frame(treatment_permanova),
          "LabNotebook/Chap4/treatment_permanova_results.csv")

# Pairwise PERMANOVA comparisons (if overall test is significant)
if(treatment_permanova$`Pr(>F)`[1] < 0.05) {
  
  cat("\n=== Running Pairwise PERMANOVA Comparisons ===\n")
  
  # Treated PMS vs Untreated PMS (both have disease_course = PPMS or SPMS)
  ms_rare_pms_only <- subset_samples(ms_rare, disease_course != "Control")
  pms_only_dist <- distance(ms_rare_pms_only, method = "wunifrac")
  pms_only_metadata <- data.frame(sample_data(ms_rare_pms_only))
  
  treated_vs_untreated <- adonis2(pms_only_dist ~ treatment_status, 
                                  data = pms_only_metadata, 
                                  permutations = 999)
  cat("\n=== Treated PMS vs Untreated PMS ===\n")
  print(treated_vs_untreated)
  
  # Healthy Control vs Untreated PMS
  ms_rare_ctrl_untreated <- subset_samples(ms_rare, 
                                           disease_course == "Control" | treatment_status == "Untreated")
  ctrl_untreated_dist <- distance(ms_rare_ctrl_untreated, method = "wunifrac")
  ctrl_untreated_metadata <- data.frame(sample_data(ms_rare_ctrl_untreated))
  
  # Create a combined variable for this comparison
  ctrl_untreated_metadata$comparison_group <- ifelse(ctrl_untreated_metadata$disease_course == "Control",
                                                     "Healthy Control",
                                                     "Untreated PMS")
  
  ctrl_vs_untreated <- adonis2(ctrl_untreated_dist ~ comparison_group, 
                               data = ctrl_untreated_metadata, 
                               permutations = 999)
  cat("\n=== Healthy Control vs Untreated PMS ===\n")
  print(ctrl_vs_untreated)
  
  # Healthy Control vs Treated PMS
  ms_rare_ctrl_treated <- subset_samples(ms_rare, 
                                         disease_course == "Control" | treatment_status == "Treated")
  ctrl_treated_dist <- distance(ms_rare_ctrl_treated, method = "wunifrac")
  ctrl_treated_metadata <- data.frame(sample_data(ms_rare_ctrl_treated))
  
  # Create a combined variable for this comparison
  ctrl_treated_metadata$comparison_group <- ifelse(ctrl_treated_metadata$disease_course == "Control",
                                                   "Healthy Control",
                                                   "Treated PMS")
  
  ctrl_vs_treated <- adonis2(ctrl_treated_dist ~ comparison_group, 
                             data = ctrl_treated_metadata, 
                             permutations = 999)
  cat("\n=== Healthy Control vs Treated PMS ===\n")
  print(ctrl_vs_treated)
  
  # Compile results
  pairwise_results <- data.frame(
    Comparison = c("Treated vs Untreated PMS", 
                   "Healthy Control vs Untreated PMS", 
                   "Healthy Control vs Treated PMS"),
    R2 = c(treated_vs_untreated$R2[1], 
           ctrl_vs_untreated$R2[1], 
           ctrl_vs_treated$R2[1]),
    F_statistic = c(treated_vs_untreated$F[1], 
                    ctrl_vs_untreated$F[1], 
                    ctrl_vs_treated$F[1]),
    p_value = c(treated_vs_untreated$`Pr(>F)`[1], 
                ctrl_vs_untreated$`Pr(>F)`[1], 
                ctrl_vs_treated$`Pr(>F)`[1])
  )
  
  # FDR correction for multiple comparisons (Benjamini-Hochberg method)
  pairwise_results$p_fdr <- p.adjust(pairwise_results$p_value, method = "BH")
  
  # Add significance indicator
  pairwise_results$significant_fdr <- ifelse(pairwise_results$p_fdr < 0.05, "Yes", "No")
  
  cat("\n=== Pairwise Comparisons Summary ===\n")
  print(pairwise_results)
  
  n_sig_unadjusted <- sum(pairwise_results$p_value < 0.05)
  n_sig_fdr <- sum(pairwise_results$p_fdr < 0.05)
  
  # Display significant comparisons
  if(n_sig_fdr > 0) {
    cat("\n=== Significant Comparisons (FDR < 0.05) ===\n")
    sig_comparisons <- pairwise_results[pairwise_results$p_fdr < 0.05, ]
    for(i in 1:nrow(sig_comparisons)) {
      cat(sig_comparisons$Comparison[i], 
          ": R² =", round(sig_comparisons$R2[i], 4),
          ", p =", round(sig_comparisons$p_value[i], 3),
          ", p_FDR =", round(sig_comparisons$p_fdr[i], 3), "\n")
    }
  } else {
    cat("\n=== No comparisons significant after FDR correction ===\n")
  }
  
} else {
  cat("\n=== Overall PERMANOVA not significant (p ≥ 0.05) ===\n")
  cat("Pairwise comparisons not performed.\n")
}

write.csv(pairwise_results, 
          "LabNotebook/Chap4/treatment_pairwise_permanova.csv", 
          row.names = FALSE)