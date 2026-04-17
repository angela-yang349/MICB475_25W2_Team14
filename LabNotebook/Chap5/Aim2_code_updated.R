# Load the packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggsignif)
library(dplyr)

# Load rarefied data
load("LabNotebook/Chap4/ms_rare_no_RRMS_ctrl.RData")


#### AIM 2: Alpha diversity across treatment types (with controls) ####
corrected_treatment_types_shannon_plot <-plot_richness(ms_rare_no_RRMS_ctrl, 
                                           x = "treatments", 
                                           measures = c("Shannon")) +
  xlab("Treatment Type") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  geom_point() +
  ggtitle("Alpha Diversity Across Treatment Types and Healthy Controls")+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
corrected_treatment_types_shannon_plot

#save the file as a .png on your local computer
#ggsave(filename = "LabNotebook/Chap5/corrected_treatment_types_shannon_plot.png",
       #corrected_treatment_types_shannon_plot,
       #height=5, width=7, dpi = 300)

## Statistical test - Kruskal-Wallis rank sum test

# Extract richness metrics and metadata
Aim2_richness_estimate <- estimate_richness(ms_rare_no_RRMS_ctrl, 
                                            measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
Aim2_alpha_samp_dat <- sample_data(ms_rare_no_RRMS_ctrl)
Aim2_alpha_samp_and_richness <- data.frame(Aim2_alpha_samp_dat, Aim2_richness_estimate)

view(Aim2_alpha_samp_and_richness)

#run Kruskal-Wallis rank sum test
kruskal_treatments <- kruskal.test( Shannon ~ treatments, data = Aim2_alpha_samp_and_richness)
kruskal_treatments

#none of these differences are significant since p=0.07782 (not less than 0.05), thus further analysis is not needed

write.csv(data.frame(
  statistic = kruskal_treatments$statistic,
  p_value = kruskal_treatments$p.value,
  method = kruskal_treatments$method,
  parameter = kruskal_treatments$parameter
), "LabNotebook/Chap5/corrected_kruskal_treatments.csv", row.names = FALSE)


#### AIM 2: Beta diversity across treatment types (with controls) ####

treatment_types_wunifrac_dist <- distance(ms_rare_no_RRMS_ctrl, method = "wunifrac")

treatment_types_wunifrac_pcoa <- ordinate(ms_rare_no_RRMS_ctrl, 
                                          method = "PCoA", 
                                          distance = treatment_types_wunifrac_dist)

# Create PCoA plot WITH healthy controls (colored by treatment type)
treatment_types_beta_plot_all <- plot_ordination(ms_rare_no_RRMS_ctrl, 
                                                 treatment_types_wunifrac_pcoa, 
                                                 color = "treatments") +
  labs(color = "Treatment Type",
       title = "Beta Diversity Across Treatment Types and Healthy Controls",
       subtitle = "PCoA with Weighted UniFrac Distance") +
  theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 8))
treatment_types_beta_plot_all

#ggsave("LabNotebook/Chap5/updated_treatment_types_wunifrac_pcoa_with_controls.png", 
       #treatment_types_beta_plot_all, 
       #width = 9, height = 6, dpi = 300)

# Create subset WITHOUT healthy controls
ms_rare_pms_treatments <- subset_samples(ms_rare_no_RRMS_ctrl, disease_course != "Control")

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

#ggsave("LabNotebook/Chap5/updated_treatment_types_wunifrac_pcoa_pms_only.png", 
       #treatment_types_beta_plot_pms, 
       #width = 9, height = 6, dpi = 300)

# PERMANOVA test (treatment type effect)
treatment_types_metadata_df <- data.frame(sample_data(ms_rare_no_RRMS_ctrl))

set.seed(123)
treatment_types_permanova <- adonis2(treatment_types_wunifrac_dist ~ treatments, 
                                     data = treatment_types_metadata_df,
                                     permutations = 999)
cat("\n=== PERMANOVA: Beta Diversity Across Treatment Types ===\n")
print(treatment_types_permanova)

# Save overall PERMANOVA results
#write.csv(as.data.frame(treatment_types_permanova),
#"LabNotebook/Chap5/updated_treatment_types_permanova_results.csv")

#differences are significant (p=0.011) so do pairwise PERMANOVA comparisons

if(treatment_types_permanova$`Pr(>F)`[1] < 0.05) {
  
  cat("\n=== Running Pairwise PERMANOVA Comparisons ===\n")
  cat("Strategy: Each treatment vs Untreated, Untreated vs Control, Each treatment vs Control\n\n")
  
  # Get list of all treatment types (excluding Untreated and Control)
  all_treatments <- unique(as.character(sample_data(ms_rare_no_RRMS_ctrl)$treatments))
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
    ms_rare_comparison <- subset_samples(ms_rare_no_RRMS_ctrl, 
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
  
  ms_rare_untreated_ctrl <- subset_samples(ms_rare_no_RRMS_ctrl, 
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
    ms_rare_comparison <- subset_samples(ms_rare_no_RRMS_ctrl, 
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
  
  # Sort by p-value
  pairwise_comparisons <- pairwise_comparisons[order(pairwise_comparisons$p_value), ]
  
  # Display summary
  cat("\n=== Pairwise PERMANOVA Summary ===\n")
  print(pairwise_comparisons, row.names = FALSE)
  
  # Count significant comparisons
  n_sig_unadjusted <- sum(pairwise_comparisons$p_value < 0.05)
  
  # Visual summary of significant comparisons (unadjusted)
  if(n_sig_unadjusted > 0) {
    cat("\n=== Significant Comparisons (p < 0.05) ===\n")
    sig_comparisons <- pairwise_comparisons[pairwise_comparisons$p_value < 0.05, ]
    for(i in 1:nrow(sig_comparisons)) {
      cat(sig_comparisons$Comparison[i], 
          ": R² =", round(sig_comparisons$R2[i], 4),
          ", p =", round(sig_comparisons$p_value[i], 4), "\n")
    }
  } else {
    cat("\n=== No comparisons were significant at p < 0.05 ===\n")
  }
  
} else {
  cat("\n=== Overall PMS-only PERMANOVA not significant (p ≥ 0.05) ===\n")
  cat("Pairwise comparisons not performed.\n")
}


#write.csv(pairwise_comparisons,
          #"LabNotebook/Chap5/updated_treatment_types_pairwise_permanova.csv",
          #row.names = FALSE)


######## AIM 2 FINAL FIGURES #########
#make sure to run everything before this to generate the figures

#relabeling of treatment names
Aim2_alpha_samp_and_richness <- Aim2_alpha_samp_and_richness %>%
  mutate(treatments = recode(treatments,
                             "control" = "Healthy Control",
                             "Dimethyl fumarate" = "Dimethyl Fumarate",
                             "Fingolimod" = "Fingolimod",
                             "Glatiramer acetate" = "Glatiramer Acetate",
                             "Interferon" = "Interferon",
                             "ocrevus(rituxan)" = "Ocrevus/Rituxan",
                             "Untreated" = "Untreated PMS"
  )) %>%
  mutate(treatments = factor(treatments, levels = c(
    "Healthy Control",
    "Dimethyl Fumarate",
    "Glatiramer Acetate",
    "Interferon",
    "Ocrevus/Rituxan",
    "Fingolimod",
    "Untreated PMS"
  )))

### Figure S2A - Alpha diversity using individual treatments
final_figS2A <- ggplot(Aim2_alpha_samp_and_richness, aes(x = treatments, y = Shannon)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill = treatments)) +
  scale_fill_manual(values = c(
    "Control" = "#e41a1c",
    "Dimethyl fumarate" = "#377eb8",
    "Fingolimod" = "#4daf4a",
    "Glatiramer acetate" = "#984ea3",
    "Interferon" = "#ff7f00",
    "Ocrevus (Rituxan)" = "#999999",
    "Untreated" = "#17334f"
  )) +
  labs(x = "Treatment", y = "Shannon Diversity Index") +
  ylim(0.0, 2.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 24),
    axis.title.y = element_text(margin = margin(r = 15)),  # adds space from axis
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20),  # extra breathing room
    legend.position = "none"
  )

final_figS2A

#New alpha plot
final_figS2A <- ggplot(Aim2_alpha_samp_and_richness, aes(x = treatments, y = Shannon)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill = treatments)) +
  labs(x = "Treatment", y = "Shannon Diversity Index") +
  ylim(0.0, 2.5) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 24),
    axis.title.y = element_text(margin = margin(r = 15)),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20),
    legend.position = "none"
  )

final_figS2A

#ggsave("LabNotebook/Chap5/updated_final_figS2A.png", final_figS2A, height = 8, width = 12)


### Figure S2B - Beta diversity using individual treatments
# Calculate percent variance
percent_var2 <- pms_treatments_wunifrac_pcoa$values$Relative_eig[1:2] * 100

axis_labels <- c(
  paste0("Axis 1 [", round(percent_var2[1], 1), "%]"),
  paste0("Axis 2 [", round(percent_var2[2], 1), "%]")
)

# Create plot with solid ellipses
final_figS2B <- plot_ordination(
  ms_rare_no_RRMS_ctrl, 
  pms_treatments_wunifrac_pcoa, 
  color = "treatments"
) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = seq(-0.5, 0.5, by = 0.25)) +
  stat_ellipse(aes(color = treatments), type = "norm", size = 0.8) +
  scale_color_manual(values = c(
    "Control" = "#e41a1c",
    "Dimethyl fumarate" = "#377eb8",
    "Fingolimod" = "#4daf4a",
    "Glatiramer acetate" = "#984ea3",
    "Interferon" = "#ff7f00",
    "Ocrevus (Rituxan)" = "#999999",
    "Untreated" = "#17334f"
  )) +
  labs(
    x = axis_labels[1],
    y = axis_labels[2], 
    color = "Treatments"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20)
  )

final_figS2B

#changing labels and order
sample_data(ms_rare_no_RRMS_ctrl)$treatments <- recode(
  sample_data(ms_rare_no_RRMS_ctrl)$treatments,
  "Control" = "Healthy Control",
  "Dimethyl fumarate" = "Dimethyl Fumarate",
  "Fingolimod" = "Fingolimod",
  "Glatiramer acetate" = "Glatiramer Acetate",
  "Interferon" = "Interferon",
  "ocrevus(rituxan)" = "Ocrevus/Rituxan",
  "Untreated" = "Untreated PMS"
)

sample_data(ms_rare_no_RRMS_ctrl)$treatments <- factor(
  sample_data(ms_rare_no_RRMS_ctrl)$treatments,
  levels = c(
    "Healthy Control",
    "Dimethyl Fumarate",
    "Glatiramer Acetate",
    "Interferon",
    "Ocrevus/Rituxan",
    "Fingolimod",
    "Untreated PMS"
  )
)

#New PcoA plot
percent_var2 <- treatment_types_wunifrac_pcoa$values$Relative_eig[1:2] * 100

axis_labels <- c(
  paste0("Axis 1 [", round(percent_var2[1], 1), "%]"),
  paste0("Axis 2 [", round(percent_var2[2], 1), "%]")
)

# Create plot with default colors
final_figS2B <- plot_ordination(
  ms_rare_no_RRMS_ctrl, 
  treatment_types_wunifrac_pcoa, 
  color = "Treatment Type"
) +
  geom_point(aes(color = treatments), size = 2) +
  scale_y_continuous(
    breaks = seq(-0.5, 0.5, by = 0.25),
    limits = c(-0.55, 0.55)
  ) +
  stat_ellipse(aes(color = treatments), type = "norm", size = 0.8) +
  labs(
    x = axis_labels[1],
    y = axis_labels[2], 
    color = "Treatment Type"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20)
  )

final_figS2B

ggsave("LabNotebook/Chap5/updated_final_figS2B.png", final_figS2B, height = 8, width = 12)
