# install package
install.packages("indicspecies")

# !/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

# set seed to account for randomness
set.seed(100)

# load data
load("ms_rare.RData")
ms_rare_no_RRMS_ctrl <- subset_samples(ms_rare, disease_course_control != "Control_RRMS")

# Indicator Species/Taxa Analysis
# glom to genus
ms_genus <- tax_glom(ms_rare_no_RRMS_ctrl, "Genus", NArm = FALSE)
ms_genus_RA <- transform_sample_counts(ms_genus, fun=function(x) x/sum(x))

# ISA
isa_ms <- multipatt(t(otu_table(ms_genus_RA)), cluster = sample_data(ms_genus_RA)$`treatments`)
summary(isa_ms)
taxtable <- tax_table(ms_rare_no_RRMS_ctrl) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of anything beyond the glomed taxa level
res <- isa_ms$sign %>% 
  rownames_to_column(var="ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value<0.05)

# view results
view(res)

#save table to github
write.csv(res, "ISA_results.csv", row.names = FALSE) 

# ===============================
# ISA grouped by treatment_status (control, treated, untreated)

isa_ms_status <- multipatt(
  t(otu_table(ms_genus_RA)), 
  cluster = sample_data(ms_genus_RA)$treatment_status
)

summary(isa_ms_status)

res_status <- isa_ms_status$sign %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value < 0.05)

# view results
View(res_status)

# save table
write.csv(res_status, "ISA_results_treatment_status.csv", row.names = FALSE)
