# Load packages
library(tidyverse)
library(phyloseq)
library(vegan)

# Load rarefied data
load("LabNotebook/Chap4/ms_rare_no_RRMS_ctrl.RData")
print(table(sample_data(ms_rare_no_RRMS_ctrl)$treatments))

sample_df <- data.frame(sample_data(ms_rare_no_RRMS_ctrl))

# ===== Group 1: healthy controls vs untreated PMS vs immunomodulators (Glatiramer acetate, Interferon (IFN), Dimethyl fumarate) vs T/B cell treatments (Ocrevus (rituxan), Fingolimod) =====
sample_df$treatment_mechanism_4grp <- case_when(
  sample_df$treatments == "Healthy Control" ~ "Healthy Control",
  sample_df$treatments == "Untreated" ~ "Untreated PMS",
  sample_df$treatments %in% c("Glatiramer acetate", "Interferon", "Dimethyl fumarate") ~ "Immunomodulators",
  sample_df$treatments %in% c("ocrevus(rituxan)", "Fingolimod") ~ "T/B Cell Therapies",
  TRUE ~ NA_character_
)

# ===== Group 2: healthy controls vs untreated PMS vs immunomodulators (Glatiramer acetate, Interferon (IFN)) vs T/B cell treatments (Ocrevus (rituxan), Fingolimod) vs DMF =====
sample_df$treatment_mechanism_5grp <- case_when(
  sample_df$treatments == "Healthy Control" ~ "Healthy Control",
  sample_df$treatments == "Untreated" ~ "Untreated PMS",
  sample_df$treatments %in% c("Glatiramer acetate", "Interferon") ~ "Immunomodulators",
  sample_df$treatments %in% c("ocrevus(rituxan)", "Fingolimod") ~ "T/B Cell Therapies",
  sample_df$treatments == "Dimethyl fumarate" ~ "DMF",
  TRUE ~ NA_character_
)

# Update phyloseq object with new groupings
sample_data(ms_rare_no_RRMS_ctrl)$treatment_mechanism_4grp <- sample_df$treatment_mechanism_4grp
sample_data(ms_rare_no_RRMS_ctrl)$treatment_mechanism_5grp <- sample_df$treatment_mechanism_5grp


### ALPHA DIVERSITY ANALYSIS ###
GROUPING <- "treatment_mechanism_4grp"

mechanism_alpha_plot <- plot_richness(ms_rare, 
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
