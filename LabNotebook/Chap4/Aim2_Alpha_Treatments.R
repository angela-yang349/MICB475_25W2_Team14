#load the packages
library(tidyverse)
library(phyloseq)
library(vegan)

#load rareified data
load("LabNotebook/Chap2/ms_rare.RData")

# Alpha diversity - Shannon Index
PMS_treatments_richness_plot <-plot_richness(ms_rare, x = "treatments", measures = c("Shannon")) +
  xlab("Treatments") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  ggtitle("Shannon Diversity of PMS Patients per Treatment and Healthy Controls")

PMS_treatments_richness_plot

#save the file as a .png on your local computer
ggsave(filename = "LabNotebook/Chap3/PMS_treatments_richness_plot.png",
       PMS_treatments_richness_plot,
       height=4, width=7)