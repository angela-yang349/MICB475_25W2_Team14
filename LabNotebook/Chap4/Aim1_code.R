#load the packages
library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)
library(picante)

# load in data
load("LabNotebook/Chap3/ms_phyloseq.RData")
load("LabNotebook/Chap3/ms_rare.RData")

# remove control samples
ms_rare_no_ctrl <- subset_samples(ms_rare, disease_course != "Control")

#### Beta diversity of SPMS vs. PPMS (no control)#####
wu_no_ctrl_distance <- distance(ms_rare_no_ctrl, method="wunifrac")

wu_no_ctrl_pcoa <- ordinate(ms_rare_no_ctrl, method="PCoA", distance=wu_no_ctrl_distance)

pms_no_ctrl_pcoa <- plot_ordination(ms_rare_no_ctrl, wu_no_ctrl_pcoa, color = "disease_course") +
  labs(col = "Disease Course") +
  theme_classic()
pms_no_ctrl_pcoa

ggsave("LabNotebook/Chap4/beta_SPMSandPPMS_pcoa.png"
       , pms_no_ctrl_pcoa
       , height=4, width=5)

PMS_data <- data.frame(sample_data(ms_rare_no_ctrl))
adonis2(wu_no_ctrl_distance ~ disease_course, data=PMS_data)

#### Alpha diversity (Shannon Index) of treated vs. untreated PMS (with control) #####
PMS_richness_plot <-plot_richness(ms_rare, x = "treatment_status", measures = c("Shannon")) +
  xlab("Treatment Status") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  ggtitle("Shannon Diversity of PMS Patients and Healthy Controls")

PMS_richness_plot

ggsave(filename = "LabNotebook/Chap4/alpha_PMS_plot.png",
       PMS_richness_plot,
       height=4, width=6)

#### Beta diversity of treated vs. untreated PMS (with control ) #####
