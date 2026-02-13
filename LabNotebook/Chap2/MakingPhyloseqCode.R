# Creating phyloseq object

# Load in the necessary packages
library(tidyverse)
library(ape)
library(phyloseq)
library(vegan)

# Load in the metadata,  OTU table, taxonomy file, and phylogenetic tree
metaFP  <- "LabNotebook/Chap1/ms_metadata_no_rrms.tsv"
meta <- read_delim(metaFP, delim="\t")

taxFP <- "LabNotebook/Chap2/taxonomy.tsv"
tax <- read_delim(taxFP, delim="\t")

treeFP <- "LabNotebook/Chap2/tree.nwk"
tree <- read.tree(treeFP)

otuFP <- "LabNotebook/Chap2/feature-table.txt"
otu <- read_delim(file=otuFP, delim="\t", skip=1)

# Adjust files to be read into a phyloseq object.

### re-format OTU table for phyloseq
# save as a matrix without the OTU IDs
otu_mat <- as.matrix(otu[,-1])
# make row names the OTU ID
rownames(otu_mat) <- otu$`#OTU ID`
# create phyloseq object
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

### re-format metadata for phyloseq
# save as dataframe without sample name
meta_df <- as.data.frame(meta[,-1])
# make row names the sample name
rownames(meta_df)<- meta$'sample-id'
# create phyloseq object
META <- sample_data(meta_df)

### re-format taxonomy for phyloseq
# save as a matrix with taxon strings separated into level of taxonomy and confidence removed
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
# remove feature IDs
tax_mat <- tax_mat[,-1]
# make row names the sample IDs
rownames(tax_mat) <- tax$`Feature ID`
# create phyloseq object
TAX <- tax_table(tax_mat)

### create phyloseq object
ms_phyloseq <- phyloseq(OTU, META, TAX, tree)

# View components of phyloseq object
# otu_table(ms_phyloseq)
# sample_data(ms_phyloseq)
# tax_table(ms_phyloseq)
# phy_tree(ms_phyloseq)

# save file
save(ms_phyloseq, file="LabNotebook/Chap2/ms_phyloseq.RData")
