library(tidyverse)

meta <- read_delim(file = "ms_metadata.tsv", delim = "\t")
manifest <- read_delim(file = "ms_manifest.tsv", delim = "\t")

#filter out samples for RRMS
filtered_meta <- subset(meta, disease_course != "RRMS")

#filter out the same samples from meta in the manifest file
filtered_manifest <- subset(manifest, `sample-id` %in% filtered_meta$`sample-id`)

#save filtered files
write.table(filtered_meta,
            "ms_metadata_no_rrms.tsv",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
write.table(filtered_manifest, 
            "ms_manifest_no_rrms.tsv", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)