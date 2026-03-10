# install package
install.packages("indicspecies")

# !/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

# load data
load("ms_phyloseq.RData")

# Indicator Species/Taxa Analysis
# glom to genus
ms_genus <- tax_glom(ms_final, "Genus", NArm = FALSE)
ms_genus_RA <- transform_sample_counts(ms_genus, fun=function(x) x/sum(x))

# ISA
isa_ms <- multipatt(t(otu_table(ms_genus_RA)), cluster = sample_data(ms_genus_RA)$`donor_status`)
summary(isa_ms)
taxtable <- tax_table(ms_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of anything beyond the glomed taxa level
res <- isa_s$sign %>% 
  rownames_to_column(var="ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value<0.05)

# view results
view(res)
