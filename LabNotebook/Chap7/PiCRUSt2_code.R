## Load libraries
#install.packages('MicrobiomeStat')
#install.packages('ggpicrust2')
#BiocManager::install("KEGGREST")
library(tidyverse)
library(phyloseq)
library(ggpicrust2)
library(stringr)

############################ data wrangling

## load in data (use non-rarified data)
load("LabNotebook/Chap3/ms_phyloseq.RData")

## remove all control
pms_phyloseq <- subset_samples(ms_phyloseq, disease_course != "Control")

## extract metadata
meta <- sample_data(pms_phyloseq) |> 
  data.frame() |> 
  rownames_to_column("sample_name")

# group columns by treatment groups
meta <- meta %>%
  mutate(group_treatment = case_when(
    treatments %in% c("Glatiramer acetate", "Interferon", "Dimethyl fumarate") ~ "Immunomodulators",
    treatments %in% c("ocrevus(rituxan)", "Fingolimod") ~ "T_B_cell",
    treatments == "Untreated" ~ "Untreated"))

## load in KO table
ko = read.delim("LabNotebook/Chap7/pred_metagenome_unstrat.tsv",
                row.names = 1,
                check.names = FALSE)

# filter for immunomodulators 
immuno_vs_unt = meta %>% filter(group_treatment %in% c("Immunomodulators", "Untreated"))
ko_immuno = ko %>% select(all_of(immuno_vs_unt$sample_name))

# filter for T and B cell
tb_vs_unt = meta %>% filter(group_treatment %in% c("T_B_cell", "Untreated"))
ko_tb = ko %>% select(all_of(tb_vs_unt$sample_name))

############################ DAA

## run DAA using LinDA for immunomodulators vs. untreated
daa_immuno = pathway_daa(abundance = ko_immuno,
                         metadata = immuno_vs_unt, 
                         group = "group_treatment", 
                         daa_method = "LinDA", 
                         select = NULL, reference = NULL)

## run DAA using LinDA for T and B cell vs. untreated
daa_tb = pathway_daa(abundance = ko_tb,
                     metadata = tb_vs_unt, 
                     group = "group_treatment", 
                     daa_method = "LinDA", 
                     select = NULL, reference = NULL)

## annotate DAA for immunomodulators
daa_immuno_annot <- pathway_annotation(pathway = "KO", 
                                       daa_results_df = daa_immuno,
                                       ko_to_kegg = TRUE)

# annotate DAA for t and b cell
daa_tb_annot <- pathway_annotation(pathway = "KO", 
                                   daa_results_df = daa_tb,
                                   ko_to_kegg = TRUE)

## save daa results
#saveRDS(daa_immuno_annot, "LabNotebook/Chap7/daa_immuno_annot.rds")
#saveRDS(daa_tb_annot, "LabNotebook/Chap7/daa_tb_annot.rds")

############################ plots! (please help me fix them :0)

glimpse(daa_immuno_annot)
glimpse(daa_tb_annot)

summary(daa_immuno_annot)
summary(daa_tb_annot)

## plot immunomodulators vs. untreated
plot_immuno <- daa_immuno_annot %>%
  filter(p_adjust < 0.05) %>%             # only significant pathways
  arrange(desc(abs(log2_fold_change))) %>%       # sort by fold change
  head(20) %>%                            # top 20 pathways
  ggplot(aes(x = reorder(str_wrap(pathway_name, width = 30), log2_fold_change), y = log2_fold_change)) +
  geom_col(fill = "steelblue") +
  coord_flip() +                          # horizontal bars
  labs(x = "Pathway", 
       y = "Log2 Fold Change", 
       title = "Immunomodulators vs Untreated") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(r = 15)),  # adds space from axis
    axis.text = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey80", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey90", size = 0.3),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20),  # extra breathing room
    legend.position = "none"
  )

plot_immuno


## plot t and b cell vs. untreated
plot_tb <- daa_tb_annot %>%
  # filter(p_adjust < 1) %>%             # only significant pathways
  arrange(desc(abs(log2_fold_change))) %>%       # sort by fold change
  head(20) %>%                            # top 20 pathways
  ggplot(aes(x = reorder(pathway_name, log2_fold_change), y = log2_fold_change)) +
  geom_col(fill = "steelblue") +
  coord_flip() +                          # horizontal bars
  labs(x = "Pathway",
       y = "Log2 Fold Change") +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(r = 15)),  # adds space from axis
    axis.text = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.margin = margin(t = 10, r = 20, b = 10, l = 20),  # extra breathing room
    legend.position = "none"
  )

plot_tb

# save plots
#ggsave("immuno_vs_untreated_picrust.png", plot = plot_immuno, width = 8, height = 6, dpi = 300)
