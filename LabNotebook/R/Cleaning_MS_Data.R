library(tidyverse)

meta <- read_delim(file = "ms_metadata.tsv", delim = "\t")
manifest <- read_delim(file = "ms_manifest.tsv", delim = "\t")

filtered_meta <- subset(meta, disease_course != "RRMS")
filtered_manifest <- subset(manifest, `sample-id` %in% filtered_meta$`sample-id`)
