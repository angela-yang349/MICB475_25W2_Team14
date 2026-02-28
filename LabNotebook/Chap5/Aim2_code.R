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

treatment_types_permanova <- adonis2(treatment_types_wunifrac_dist ~ treatments, 
                                     data = treatment_types_metadata_df,
                                     permutations = 999)
cat("\n=== PERMANOVA: Beta Diversity Across Treatment Types ===\n")
print(treatment_types_permanova)

#differences are significant (p=0.002) so do pairwise PERMANOVA comparisons


