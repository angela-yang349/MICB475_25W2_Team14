# February 10th, 2026
## Agenda
  * Completed work:
    * Filtered metadata and manifest files to remove RRMS patients
    * Ran QIIME2 pipeline on filtered data and generated table.qzv and rep-seq.qzv
  * Upcoming work:
    * Diversity metrics
      * Keep SPMS and PPMS separate then look at different treatments within each group vs. non-treated
      * Could bin treatments that have the same significance
    * Team proposal 
  * Questions:
    * How should we determine sampling depth?
    * Team proposal

## Discussion
1. Proposal Tips From Ritu
  * Title of proposal should start with -ing word
  * Mention previous papers
  * Talk about the different types of MS (RRMS, PPMS)
  * Research question section, put question in bold and mention hypothesis
    * Hypothesis should be based on previous papers and in a formal tone
  * Aims should be titled like a sentence (not just “diversity analysis”)
  * Mention which statistic tests we are using in the table for the proposal
  * Key: look at example proposals as a guideline
2. Current Progress
  * In process of picking sampling depth → 55
  * table.qzv and rep-seq.qzv are in the repository
  * Quality of sequences is good, no trimming needed
  * Working on rarefaction part
3. Meeting During Reading Break?
  * Either a Zoom meeting later in the week or we email Ritu questions from our own separate meeting
4. Future Analysis
  * Indicator species - looking for species that are present in one group and not present in another group, generally no threshold, look at how the abundance of this species changes, down/upregulation of species
  * Can still move forward with rest of analysis if diversity metrics show no significance
  * QIIME processing might not need to be part of proposal (examples have done aim 1 = diversity, 1a is QIIME)
  * Metadata - sample IDs, characteristics like age, which disease they have
  * Manifest -  sequences, corresponds to an ID
  * Remove mitochondria and chloroplasts then do phylogenetic tree and then rarefaction done in QIIME

## Conclusions
  * Continue working on proposal, basic data processing in QIIME2, and getting organized on GitHub

## To-do
- [ ] Create files needed for phyloseq object: metadata, rooted tree, filtered table, taxonomy
- [ ] Organize our GitHub and lab notebook section
- [ ] Work on proposal


