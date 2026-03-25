# Load libraries
library(tidyverse)
library(phyloseq)
library(ggpicrust2)

# load in data (use non-rarified data)
load("LabNotebook/Chap3/ms_phyloseq.RData")

# remove all control
pms_phyloseq <- subset_samples(ms_phyloseq, disease_course != "Control")

