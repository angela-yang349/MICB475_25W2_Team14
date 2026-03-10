# install package
install.packages("indicspecies")

# !/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

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

# save as a bar chart to github
library(ggplot2)

ggplot(res, aes(x = reorder(OTU, stat), y = stat)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() -> p

ggsave("ISA_barplot.png", p, width = 10, height = 7)

