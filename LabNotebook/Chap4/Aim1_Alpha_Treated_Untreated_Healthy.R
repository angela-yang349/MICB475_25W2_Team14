#load the packages
library(tidyverse)
library(phyloseq)
library(vegan)

#load rareified data
load("LabNotebook/Chap2/ms_rare.RData")

# Alpha diversity - Shannon Index
PMS_richness_plot <-plot_richness(ms_rare, x = "treatment_status", measures = c("Shannon")) +
  xlab("Treatment Status") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  ggtitle("Shannon Diversity of PMS Patients and Healthy Controls")

PMS_richness_plot

#save the file as a .png on your local computer
ggsave(filename = "LabNotebook/Chap3/PMS_richness_plot.png",
       PMS_richness_plot,
       height=4, width=6)