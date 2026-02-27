#load the packages
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggsignif)

#load rareified data
load("LabNotebook/Chap2/ms_rare.RData")

# Alpha diversity - Shannon Index
PMS_treatments_richness_plot <-plot_richness(ms_rare, x = "treatments", measures = c("Shannon")) +
  xlab("Treatments") +
  ylab("Shannon Diversity Index") +
  geom_boxplot() +
  geom_point() +
  ggtitle("Shannon Diversity of PMS Patients per Treatment and Healthy Controls")+
  theme_classic()

PMS_treatments_richness_plot

#save the file as a .png on your local computer
#ggsave(filename = "LabNotebook/Chap5/PMS_treatments_richness_plot.png",
       #PMS_treatments_richness_plot,
       #height=4, width=7)

### Statistical test - Kruskal-Wallis rank sum test

#first extract info
Aim2_richness_estimate <- estimate_richness(ms_rare, 
                                            measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
Aim2_alpha_samp_dat <- sample_data(ms_rare)
Aim2_alpha_samp_and_richness <- data.frame(Aim2_alpha_samp_dat, Aim2_richness_estimate)

view(Aim2_alpha_samp_and_richness)

#run Kruskal-Wallis rank sum test
kruskal_treatments <- kruskal.test( Shannon ~ treatments, data = Aim2_alpha_samp_and_richness)
kruskal_treatments

#none of these differences are significant since p=0.18 (not less than 0.05), thus further analysis is not needed
