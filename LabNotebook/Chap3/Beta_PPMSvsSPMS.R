library(phyloseq)
library(ape)
library(tidyverse)
library(picante)

# load in data
load("LabNotebook/Chap2/ms_phyloseq.RData")
load("LabNotebook/Chap2/ms_rare.RData")

# remove control samples
ms_rare_no_ctrl <- subset_samples(ms_rare, disease_course != "Control")

#### Beta diversity #####
wu_no_ctrl_distance <- distance(ms_rare_no_ctrl, method="wunifrac")

wu_no_ctrl_pcoa <- ordinate(ms_rare_no_ctrl, method="PCoA", distance=wu_no_ctrl_distance)

pms_no_ctrl_pcoa <- plot_ordination(ms_rare_no_ctrl, wu_no_ctrl_pcoa, color = "disease_course") +
  labs(col = "Disease Course") +
  theme_classic()
pms_no_ctrl_pcoa

ggsave("LabNotebook/Chap3/pms_no_control_pcoa.png"
       , pms_no_ctrl_pcoa
       , height=4, width=5)

PMS_data <- data.frame(sample_data(ms_rare_no_ctrl))
adonis2(wu_no_ctrl_distance ~ disease_course, data=PMS_data)
