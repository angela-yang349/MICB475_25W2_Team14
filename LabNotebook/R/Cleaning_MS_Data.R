library(tidyverse)

meta <- read_delim(file = "ms_metadata.tsv", delim = "\t")
filtered_meta <- subset(meta, disease_course != "RRMS")
