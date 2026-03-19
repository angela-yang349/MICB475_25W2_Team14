# Group 2
# Compare grouped treatments: GA/IFN/DMF vs Ocrevus/Fingolimod

library(tidyverse)
library(phyloseq)
library(indicspecies)
install.packages("kableExtra")
library(kableExtra)

set.seed(100)

# load data
load("ms_rare.RData")
ms_rare_no_RRMS_ctrl <- subset_samples(ms_rare, disease_course_control != "Control_RRMS")

# glom to genus
ms_genus <- tax_glom(ms_rare_no_RRMS_ctrl, "Genus", NArm = FALSE)
ms_genus_RA <- transform_sample_counts(ms_genus, fun=function(x) x/sum(x))

# ===============================
# GROUPING: 2 groups based on mechanism/effect
# ===============================

sample_data(ms_genus_RA)$grouped_by_mechanism <- case_when(
  sample_data(ms_genus_RA)$treatments %in% c(
    "Glatiramer acetate",
    "Interferon",
    "Dimethyl fumarate"
  ) ~ "Immunomodulators",
  
  sample_data(ms_genus_RA)$treatments %in% c(
    "ocrevus(rituxan)",
    "Fingolimod"
  ) ~ "Lymphocyte",
  
  TRUE ~ NA_character_
)

# subset samples without NA
ms_group2 <- subset_samples(ms_genus_RA, !is.na(grouped_by_mechanism))

# run ISA
isa_group2 <- multipatt(
  t(otu_table(ms_group2)), 
  cluster = sample_data(ms_group2)$grouped_by_mechanism
)

summary(isa_group2)

taxtable <- tax_table(ms_rare_no_RRMS_ctrl) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ASV")

res_group2 <- isa_group2$sign %>% 
  rownames_to_column(var = "ASV") %>% 
  left_join(taxtable) %>%
  filter(p.value < 0.05)

View(res_group2)
write.csv(res_group2, "ISA_results_group2.csv", row.names = FALSE)

html_group2 <- res_group2 %>%
  kable(format = "html", table.attr = "border='1' style='border-collapse: collapse;'", escape = TRUE) %>%
  kable_styling(full_width = FALSE, position = "left")

save_kable(html_group2, "ISA_results_group2.html")