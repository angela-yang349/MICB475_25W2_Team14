# Load packages
library(tidyverse)
library(phyloseq)
library(vegan)

# Load rarefied data
load("LabNotebook/Chap4/ms_rare_no_RRMS_ctrl.RData")
print(table(sample_data(ms_rare_no_RRMS_ctrl)$treatments))

sample_df <- data.frame(sample_data(ms_rare_no_RRMS_ctrl))

# Group: healthy controls vs untreated PMS vs immunomodulators (Glatiramer acetate, Interferon (IFN), Dimethyl fumarate) vs T/B cell treatments (Ocrevus (rituxan), Fingolimod) =====
sample_df$treatment_mechanism_4grp <- case_when(
  sample_df$treatments == "Control" ~ "Healthy Control",
  sample_df$treatments == "Untreated" ~ "Untreated PMS",
  sample_df$treatments %in% c("Glatiramer acetate", "Interferon", "Dimethyl fumarate") ~ "Immunomodulators",
  sample_df$treatments %in% c("ocrevus(rituxan)", "Fingolimod") ~ "T/B Cell Therapies",
  TRUE ~ NA_character_
)

# Update phyloseq object with new groupings
sample_data(ms_rare_no_RRMS_ctrl)$treatment_mechanism_4grp <- sample_df$treatment_mechanism_4grp


### ALPHA DIVERSITY ANALYSIS ###
GROUPING <- "treatment_mechanism_4grp"

mechanism_alpha_plot <- plot_richness(ms_rare_no_RRMS_ctrl, 
                                      x = GROUPING, 
                                      measures = c("Shannon")) +
  xlab("Treatment Mechanism") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  geom_point() +
  ggtitle("Alpha Diversity Across Treatment Mechanism Groups") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
mechanism_alpha_plot

ggsave(paste0("LabNotebook/Chap5/grouped_treatments_shannon_plot.png"),
       mechanism_alpha_plot,
       height = 5, width = 8)

## Statistical test - Kruskal-Wallis rank sum test

# Extract richness metrics and metadata
mechanism_richness <- estimate_richness(ms_rare_no_RRMS_ctrl, measures = "Shannon")
mechanism_metadata <- data.frame(sample_data(ms_rare_no_RRMS_ctrl))
mechanism_alpha_data <- data.frame(mechanism_metadata, mechanism_richness)

#run Kruskal-Wallis rank sum test
mechanism_alpha_kruskal <- kruskal.test(as.formula(paste("Shannon ~", GROUPING)), 
                                        data = mechanism_alpha_data)
mechanism_alpha_kruskal

write.csv(data.frame(
  statistic = mechanism_alpha_kruskal$statistic,
  p_value = mechanism_alpha_kruskal$p.value,
  method = mechanism_alpha_kruskal$method,
  parameter = mechanism_alpha_kruskal$parameter
), paste0("LabNotebook/Chap5/grouoped_treatments_alpha_kruskal_", GROUPING, ".csv"), 
row.names = FALSE)


### BETA DIVERSITY ANALYSIS ###

# Calculate weighted UniFrac distance
mechanism_wunifrac_dist <- distance(ms_rare_no_RRMS_ctrl, method = "wunifrac")

# PCoA ordination
mechanism_wunifrac_pcoa <- ordinate(ms_rare_no_RRMS_ctrl, 
                                    method = "PCoA", 
                                    distance = mechanism_wunifrac_dist)

# Create PCoA plot WITH healthy controls
mechanism_beta_plot <- plot_ordination(ms_rare_no_RRMS_ctrl, 
                                           mechanism_wunifrac_pcoa, 
                                           color = GROUPING) +
  labs(color = "Treatment Mechanism",
       title = "Beta Diversity: Treatment Mechanism Groups",
       subtitle = "PCoA with Weighted UniFrac Distance") +
  theme_classic() +
  theme(legend.position = "right")
mechanism_beta_plot

ggsave(paste0("LabNotebook/Chap5/grouped_treatments_beta_pcoa.png"), 
       mechanism_beta_plot, 
       width = 10, 
       height = 7, 
       dpi = 300)

# PERMANOVA test
mechanism_metadata_df <- data.frame(sample_data(ms_rare_no_RRMS_ctrl))

mechanism_permanova <- adonis2(as.formula(paste("mechanism_wunifrac_dist ~", GROUPING)), 
                               data = mechanism_metadata_df,
                               permutations = 999)
mechanism_permanova

## PAIRWISE PERMANOVA (IF SIGNIFICANT) ##
if(pms_mechanism_permanova$`Pr(>F)`[1] < 0.05) {
  
  cat("\n===== PAIRWISE PERMANOVA: TREATMENT MECHANISMS =====\n")
  
  # Get unique mechanism groups (excluding NA)
  mechanism_groups <- unique(sample_data(ms_rare_no_RRMS_ctrl)[[GROUPING]])
  mechanism_groups <- mechanism_groups[!is.na(mechanism_groups)]
  
  cat("Groups:", paste(mechanism_groups, collapse = ", "), "\n\n")
  
  # Initialize results
  pairwise_results <- data.frame(
    Comparison = character(),
    R2 = numeric(),
    F_statistic = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Perform all pairwise comparisons
  for(i in 1:(length(mechanism_groups)-1)) {
    for(j in (i+1):length(mechanism_groups)) {
      
      group1 <- mechanism_groups[i]
      group2 <- mechanism_groups[j]
      
      # Subset to these two groups
      ms_pair <- subset_samples(ms_rare_no_RRMS_ctrl, 
                                get(GROUPING) %in% c(group1, group2))
      
      # Check sample size
      if(nsamples(ms_pair) < 3) {
        cat(group1, "vs", group2, ": Insufficient samples, skipping\n")
        next
      }
      
      # Calculate distance
      pair_dist <- distance(ms_pair, method = "wunifrac")
      pair_metadata <- data.frame(sample_data(ms_pair))
      
      # Run PERMANOVA
      pair_result <- adonis2(as.formula(paste("pair_dist ~", GROUPING)),
                             data = pair_metadata,
                             permutations = 999)
      
      # Store results
      pairwise_results <- rbind(pairwise_results,
                                data.frame(
                                  Comparison = paste(group1, "vs", group2),
                                  R2 = pair_result$R2[1],
                                  F_statistic = pair_result$F[1],
                                  p_value = pair_result$`Pr(>F)`[1]
                                ))
      
      cat(group1, "vs", group2, ": R² =", round(pair_result$R2[1], 4),
          ", p =", pair_result$`Pr(>F)`[1], "\n")
    }
  }
  
  # Sort by p-value
  pairwise_results <- pairwise_results[order(pairwise_results$p_value), ]
  
  cat("\n=== Pairwise PERMANOVA Summary ===\n")
  print(pairwise_results, row.names = FALSE)
  
  # Count significant
  n_sig <- sum(pairwise_results$p_fdr < 0.05)
  cat("\nSignificant comparisons (FDR < 0.05):", n_sig, "of", nrow(pairwise_results), "\n")
  
} else {
  cat("\n=== Overall PERMANOVA not significant, skipping pairwise comparisons ===\n")
}

  