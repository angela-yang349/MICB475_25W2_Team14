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

# Remove RRMS control samples
ms_rare_no_RRMS_ctrl <- subset_samples(ms_rare, disease_course_control != "Control_RRMS")
#save(ms_rare_no_RRMS_ctrl, file= "LabNotebook/Chap4/ms_rare_no_RRMS_ctrl.RData") 

#### AIM 1 PRELIMINARY: Beta diversity SPMS vs PPMS (no controls) ####

# Remove control samples for PMS-only analyses
ms_rare_pms_only <- subset_samples(ms_rare_no_RRMS_ctrl, disease_course != "Control")

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
corrected_treatment_alpha_plot <- plot_richness(ms_rare_no_RRMS_ctrl, 
                                      x = "treatment_status", 
                                      measures = c("Shannon")) +
  xlab("Treatment Status") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  geom_point() +
  ggtitle("Alpha Diversity: Treated vs Untreated PMS and Controls") +
  theme_classic()
corrected_treatment_alpha_plot

#ggsave(filename = "LabNotebook/Chap4/corrected_treatment_shannon_boxplot.png",
       #corrected_treatment_alpha_plot,
       #height = 4, width = 6)

### Statistical test - Kruskal-Wallis rank sum test
#first extract info
Aim1_richness_estimate <- estimate_richness(ms_rare_no_RRMS_ctrl, 
                                            measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
Aim1_alpha_samp_dat <- sample_data(ms_rare_no_RRMS_ctrl)
Aim1_alpha_samp_and_richness <- data.frame(Aim1_alpha_samp_dat, Aim1_richness_estimate)

view(Aim1_alpha_samp_and_richness)

#run Kruskal-Wallis rank sum test
kruskal_treatment_status <- kruskal.test( Shannon ~ treatment_status, data = Aim1_alpha_samp_and_richness)
kruskal_treatment_status

#none of these differences are significant since p=0.06375 (not less than 0.05), thus further analysis is not needed



#### AIM 1: Beta diversity treated vs untreated PMS (with controls) ####
treatment_wunifrac_dist <- distance(ms_rare_no_RRMS_ctrl, method = "wunifrac")

treatment_wunifrac_pcoa <- ordinate(ms_rare_no_RRMS_ctrl, 
                                    method = "PCoA", 
                                    distance = treatment_wunifrac_dist)

treatment_beta_plot <- plot_ordination(ms_rare_no_RRMS_ctrl, 
                                       treatment_wunifrac_pcoa, 
                                       color = "treatment_status") +
  labs(color = "Treatment Status", 
       title = "Beta Diversity: Treated vs Untreated PMS vs. Healthy Controls",
       subtitle = "PCoA with Weighted UniFrac Distance") +
  theme_classic() +
  theme(legend.position = "right")
treatment_beta_plot

#ggsave("LabNotebook/Chap4/updated_treatment_wunifrac_pcoa.png", 
       #treatment_beta_plot, width = 8, height = 6, dpi = 300)

# PERMANOVA test (treatment effect)
metadata_all_df <- data.frame(sample_data(ms_rare_no_RRMS_ctrl))

set.seed(123)
treatment_permanova <- adonis2(treatment_wunifrac_dist ~ treatment_status, 
                               data = metadata_all_df,
                               permutations = 999)
print(treatment_permanova)

# Save overall PERMANOVA results
#write.csv(as.data.frame(treatment_permanova),
          #"LabNotebook/Chap4/updated_treatment_permanova_results.csv")

# Pairwise PERMANOVA comparisons (if overall test is significant)
if(treatment_permanova$`Pr(>F)`[1] < 0.05) {
  
  cat("\n=== Running Pairwise PERMANOVA Comparisons ===\n")
  
  # Treated PMS vs Untreated PMS (both have disease_course = PPMS or SPMS)
  ms_rare_pms_only <- subset_samples(ms_rare_no_RRMS_ctrl, disease_course != "Control")
  pms_only_dist <- distance(ms_rare_pms_only, method = "wunifrac")
  pms_only_metadata <- data.frame(sample_data(ms_rare_pms_only))
  
  treated_vs_untreated <- adonis2(pms_only_dist ~ treatment_status, 
                                  data = pms_only_metadata, 
                                  permutations = 999)
  cat("\n=== Treated PMS vs Untreated PMS ===\n")
  print(treated_vs_untreated)
  
  # Healthy Control vs Untreated PMS
  ms_rare_ctrl_untreated <- subset_samples(ms_rare_no_RRMS_ctrl, 
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
  ms_rare_ctrl_treated <- subset_samples(ms_rare_no_RRMS_ctrl, 
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
  
  
} else {
  cat("\n=== Overall PERMANOVA not significant (p ≥ 0.05) ===\n")
  cat("Pairwise comparisons not performed.\n")
}

#write.csv(pairwise_results, 
          #"LabNotebook/Chap4/updated_treatment_pairwise_permanova.csv", 
          #row.names = FALSE)

######## AIM 1 FINAL FIGURES #########
#make sure to run everything before this to generate the figures

### Figure 1A - Alpha diversity using treatment status

# Calculate alpha diversity 
alphadiv <- estimate_richness(ms_rare_no_RRMS_ctrl, measures = c("Shannon"))

# Merge with metadata
samp_dat <- sample_data(ms_rare_no_RRMS_ctrl)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

samp_dat_wdiv$treatment_status <- factor(
  samp_dat_wdiv$treatment_status,
  levels = c("Control", "Treated", "Untreated"),
  labels = c("Healthy Control", "Treated PMS", "Untreated PMS")
)

# Alpha diversity plot (no p-value annotations)
final_fig1A <- ggplot(samp_dat_wdiv, aes(x = treatment_status, y = Shannon)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(aes(fill = treatment_status)) +
  scale_fill_manual(values = c(
    "Healthy Control" = "#e31a1c",
    "Treated PMS" = "#33a02c",
    "Untreated PMS" =  "#1f78b4"
  )) +
  labs(x = "Treatment Status", y = "Shannon Diversity Index") +
  ylim(0.4, 2.1) +
  theme_classic(base_size = 16) +  # increase overall base font
  theme(
    axis.title = element_text(size = 14),   # axis labels
    axis.text = element_text(size = 11),    # tick labels
    legend.position = "none"
  )

final_fig1A

# ggsave("LabNotebook/Chap4/new_final_fig1A.png", final_fig1A, height = 8, width = 14)


### Figure 1B - Beta diversity using treatment status
# Extract percent variance explained
percent_var <- treatment_wunifrac_pcoa$values$Relative_eig[1:2] * 100
axis_labels <- c(
  paste0("Axis 1 [", round(percent_var[1], 1), "%]"),
  paste0("Axis 2 [", round(percent_var[2], 1), "%]")
)

# Change x axis category names
sample_data(ms_rare_no_RRMS_ctrl)$treatment_status <- factor(
  sample_data(ms_rare_no_RRMS_ctrl)$treatment_status,
  levels = c("Control", "Treated", "Untreated"),
  labels = c("Healthy Control", "Treated PMS", "Untreated PMS")
)

# Create plot with solid ellipses
final_fig1B <- plot_ordination(ms_rare_no_RRMS_ctrl, 
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

final_fig1B

# ggsave("LabNotebook/Chap4/new_final_fig1B.png", final_fig1B, height = 8, width = 14)
