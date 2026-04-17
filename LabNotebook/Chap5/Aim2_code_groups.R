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

GROUPING <- "treatment_mechanism_4grp"

### ALPHA DIVERSITY ANALYSIS ###
mechanism_alpha_plot <- plot_richness(ms_rare_no_RRMS_ctrl, 
                                      x = GROUPING, 
                                      measures = c("Shannon")) +
  aes(fill = .data[[GROUPING]]) +  
  xlab("Treatment Mechanism") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Healthy Control" = "#e31a1c",
    "Untreated PMS" =  "#1f78b4",
    "Immunomodulators" = "#800080",
    "T/B Cell Therapies" = "#8A9A5B"
  )) +
  labs(fill = "Treatment Group") +
  geom_point() +
  ggtitle("Alpha Diversity Across Treatment Mechanism Groups") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),   # ↓ make y-axis numbers smaller
    axis.title = element_text(size = 10)
  )
mechanism_alpha_plot

#new plot
mechanism_alpha_plot <- plot_richness(ms_rare_no_RRMS_ctrl, 
                                      x = GROUPING, 
                                      measures = c("Shannon")) +
  aes(fill = .data[[GROUPING]]) +  
  xlab("Treatment Mechanism") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  labs(fill = "Treatment Group") +
  geom_point() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

mechanism_alpha_plot

ggsave(paste0("LabNotebook/Chap5/new_grouped_treatments_shannon_plot.png"),
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

#pairwise comparisons
pairwise_results <- pairwise.wilcox.test(
  mechanism_alpha_data$Shannon,
  mechanism_alpha_data$treatment_mechanism_4grp
)

pairwise_results

#no significant results

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

# Get unique mechanism groups (excluding NA)
mechanism_groups <- unique(sample_data(ms_rare_no_RRMS_ctrl)[[GROUPING]])
mechanism_groups <- mechanism_groups[!is.na(mechanism_groups)]

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
    
    cat("Comparing:", group1, "vs", group2, "\n")
    
    # Subset to these two groups
    ms_pair <- subset_samples(ms_rare_no_RRMS_ctrl, 
                              treatment_mechanism_4grp %in% c(group1, group2))
    
    # Check sample size
    n_samples <- nsamples(ms_pair)
    if(n_samples < 3) {
      cat("  Insufficient samples (n =", n_samples, "), skipping\n\n")
      next
    }
    
    # Calculate distance for this pair
    pair_dist <- distance(ms_pair, method = "wunifrac")
    pair_metadata <- data.frame(sample_data(ms_pair))
    
    # Run PERMANOVA for this pair
    pair_result <- adonis2(pair_dist ~ treatment_mechanism_4grp,
                           data = pair_metadata,
                           permutations = 999)
    
    # Extract and store results
    pairwise_results <- rbind(pairwise_results,
                              data.frame(
                                Comparison = paste(group1, "vs", group2),
                                R2 = pair_result$R2[1],
                                F_statistic = pair_result$F[1],
                                p_value = pair_result$`Pr(>F)`[1]
                              ))
    
    cat("  R² =", round(pair_result$R2[1], 4), 
        ", F =", round(pair_result$F[1], 2),
        ", p =", pair_result$`Pr(>F)`[1], "\n\n")
  }
}

# Add significance indicator (no FDR correction)
pairwise_results$significant <- ifelse(pairwise_results$p_value < 0.05, "Yes", "No")

# Sort by p-value
pairwise_results <- pairwise_results[order(pairwise_results$p_value), ]

# Display results
cat("\n=== PAIRWISE PERMANOVA SUMMARY ===\n")
print(pairwise_results, row.names = FALSE)

write.csv(pairwise_results,
          "LabNotebook/Chap5/grouped_treatments_pairwise_permanova.csv",
          row.names = FALSE)

# Recreate PCoA plot with stat_ellipse to show groupings
mechanism_beta_plot_ellipses <- plot_ordination(ms_rare_no_RRMS_ctrl, 
                                                mechanism_wunifrac_pcoa, 
                                                color = GROUPING) +
  stat_ellipse(aes(group = treatment_mechanism_4grp), 
               type = "norm", 
               linetype = 2, 
               size = 0.8,
               alpha = 0.3) +
  labs(color = "Treatment Mechanism",
       title = "Beta Diversity: Treatment Mechanism Groups",
       subtitle = "PCoA with Weighted UniFrac Distance\nEllipses show 95% confidence intervals",
       caption = "Significant pairwise differences (p < 0.05):\nImmunomodulators vs Control (p=0.022)\nT/B Cell vs Immunomodulators (p=0.038)") +
  theme_classic() +
  theme(legend.position = "right",
        plot.caption = element_text(hjust = 0, size = 9))
mechanism_beta_plot_ellipses

ggsave("LabNotebook/Chap5/grouped_treatments_beta_pcoa_with_ellipses.png", 
       mechanism_beta_plot_ellipses, 
       width = 11, 
       height = 8, 
       dpi = 300)

# Final Figure PCoA plot

# Step 1: extract coordinates and create axis labels
ordination_df <- plot_ordination(ms_rare_no_RRMS_ctrl, mechanism_wunifrac_pcoa, type = "samples", justDF = TRUE)

axis_labels <- c(
  paste0("Axis 1 [", round(percent_var[1], 1), "%]"),
  paste0("Axis 2 [", round(percent_var[2], 1), "%]")
)

# Step 2: create plot from scratch
mechanism_beta_plot_ellipses <- ggplot(ordination_df, aes(x = Axis.1, y = Axis.2, color = .data[[GROUPING]])) +
  # Draw ellipses
  stat_ellipse(aes(group = .data[[GROUPING]]),
               type = "norm",
               linetype = "solid",
               size = 0.8,
               alpha = 0.6) +
  # Draw points manually (you can now set size super small)
  geom_point(size = 1, alpha = 0.7) +
  # Labels
  labs( x = axis_labels[1],
        y = axis_labels[2], 
        color = "Treatment Group",
        title = "Beta Diversity: Treatment Mechanism Groups",
        subtitle = "PCoA with Weighted UniFrac Distance\nEllipses show 95% confidence intervals",
        caption = "Significant pairwise differences (p < 0.05):\nImmunomodulators vs Control (p=0.022)\nT/B Cell vs Immunomodulators (p=0.038)") +
  # Custom colors
  scale_color_manual(values = c(
    "Healthy Control" = "#e31a1c",
    "Untreated PMS" = "#1f78b4",
    "Immunomodulators" = "#800080",
    "T/B Cell Therapies" = "#8A9A5B"
  )) +
  # Expand axes so points are less clustered
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  # Fixed aspect ratio
  coord_fixed(ratio = 1) +
  # Clean theme
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.caption = element_text(hjust = 0, size = 9),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

mechanism_beta_plot_ellipses

# Final Figure 2C
mechanism_alpha_plot <- plot_richness(ms_rare_no_RRMS_ctrl, 
                                      x = GROUPING, 
                                      measures = c("Shannon")) +
  aes(fill = .data[[GROUPING]]) +  
  xlab("Treatment Mechanism") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  scale_fill_manual(values = c(
    "Healthy Control" = "#e31a1c",
    "Untreated PMS" =  "#1f78b4",
    "Immunomodulators" = "#800080",
    "T/B Cell Therapies" = "#8A9A5B"
  )) +
  labs(fill = "Treatment Group") +
  geom_point() +
  ggtitle("Alpha Diversity Across Treatment Mechanism Groups") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),   # ↓ make y-axis numbers smaller
    axis.title = element_text(size = 10)
  )
mechanism_alpha_plot

# Final Figure 2D

# Extract percent variance explained for axis labels
percent_var <- treatment_wunifrac_pcoa$values$Relative_eig[1:2] * 100
axis_labels <- c(
  paste0("Axis 1 [", round(percent_var[1], 1), "%]"),
  paste0("Axis 2 [", round(percent_var[2], 1), "%]")
)

sample_data(ms_rare_no_RRMS_ctrl)$treatment_status <- factor(
  sample_data(ms_rare_no_RRMS_ctrl)$treatment_status,
  levels = c("Control", "Treated", "Untreated"),
  labels = c("Healthy Control", "Treated PMS", "Untreated PMS")
)

# Create plot with solid ellipses
treatment_beta_plot <- plot_ordination(ms_rare_no_RRMS_ctrl, 
                                       treatment_wunifrac_pcoa, 
                                       color = "treatment_status") +
  geom_point(size = 2) +
  stat_ellipse(type = "norm", size = 0.8) +   # solid ellipse lines
  scale_color_manual(values = c(
    "Healthy Control" = "#e31a1c",
    "Treated PMS" = "#33a02c",
    "Untreated PMS" = "#1f78b4"
  )) +
  labs(
    x = axis_labels[1],
    y = axis_labels[2], 
    color = "Treatment Status"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

treatment_beta_plot
#export plot using preview
