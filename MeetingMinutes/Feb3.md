# February 3rd, 2026
## Agenda
* Discussion of project topics:
  * Vaginosis dataset
    * Very small dataset
    * May not be amplicon sequencing data because it says WGA
  * IVF dataset
    * New dataset so may have difficulty in identifying samples and metadata categories are not clear about which samples came from which patient
    * Sample collection is not consistent across patients
    * Could have potential contamination
  * MS dataset (already exists in server)
    * Dataset is massive (will take a long time to process)
    * Will need to recreate manifest and metadata without the rrms patients
    * Sample size concern (might be too small) with certain treatments

## Conclusions
* Collectively decided on MS dataset, looking into treatment differences in severe MS patients (SPMS and PPMS)

## Discussion
#### Ms Project Pipeline:
1. Processing
    * Remove RRMS samples from metadata and manifest files
2. QIIME pipeline
3. Diversity metrics
    * Keep SPMS and PPMS separate then look at different treatments within each group vs. non-treated
    * Could bin treatments that have the same significance
4. Core microbiome
    * gives venn diagram of microbes that are different and shared between groups
5. Indicator taxa
    * gives table of which microbes are strongly associated with treatment vs. non-treatment in both groups
    * can do literature search to see if the microbes are pathogenic or beneficial
6. Deseq (depends if shannon’s diversity or evenness is significant) or functional analysis (only if there is significant indicator taxa of microbes that are very functionally different and helps to predict which pathways are up- or downregulated)

## To-do
1. For data processing:
- [ ] In R, before importing the data, first sort the MS dataset from the server and copy it into our folder (keep it in .tsv file format)
- [ ] Then match the sorted sample IDs to filter for the same samples in the manifest file
2. For project 2 proposal:
- [ ] Make sure to mention and describe the different stages of MS and why we are only focusing on the 2
- [ ] Talk about what previous literature has shown on significance of treatment results
