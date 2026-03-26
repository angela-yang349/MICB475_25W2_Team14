# March 24th, 2026
## Agenda
* Completed work:
  * Redid alpha and beta diversity for Aim 1 and 2 after removal of RRMS household controls
  * Core microbiome venn diagrams
  * Indicator species analyses
  * DESeq2 Analysis - Volcano plots
  * Made improvements to figures for Aim 1 and Aim 2
  
* Upcoming work:
  * Finalizing our figures

* Questions:
  * For DESeq2 analysis, only able to find significant ASVs if I use the raw phyloseq object before rarefaction. If I use the rarefied data, I get not significance at all and am unable to run the DESeq2 analysis
  * Should we re-run our Aim2 analysis (beta and alpha) again but using the same grouped treatments from the core microbiome and ISA?
  * Should our preliminary test be a supplementary figure?
  * How do we edit the Venn diagram to make it ready for the paper?
  * What grouping should we go with for the rest of our analyses?
    * Not enough data for Ocrevus (Rituxan) - only 2 points

## Conclusions
* Only alpha and beta diversity use rarified data as an input
* Overall story = nuanced changes are not translating to community changes. Can we talk about the ones that are being impacted
* Consider the treatments that give different results and whether past papers discuss it
* ISA Table → discard (will not be making a figure for this)
    * Only keep stat value >0.7 and are unique to one condition
    * So we dont have any 😔

## Discussion
* Finalizing figures:
  * Figure 1 - Alpha and Beta diversity
    * Get rid of “treatment status” on x-axis for panel A
    * Solid lines for panel B
  * Include the new grouped treatments from Aim2 (alpha and beta) as panel C and D
    * Add significance (astrix) to panel C alpha diversity
    * Change colour so that it is more distinguishable (Make rituxan colour #ffea00)
    * Redo Aim 2 with groupings
    * Mention how sample size for some treatments were very low 
* Figure 2 - Core Microbiome
  * Evelyn prefers 2 grouped Venn diagram from core microbiome
    * Immunodulators (15), T and B cell (16)
    * Investigate into taxa that are only part of T and B cell (4) and immunomodulator (1) and also look into that taxa that is shared between immunomodulator and T and B cells (1)
    * Make font bigger
    * Go back and match to taxa table
* Figure 3 - DESeq2 Volcano plot
  * Annotate/mention the specific taxa that show “large effect” → will be talked about more in the discussion
  * Plot treatment groups vs untreated (individually)
    * Immunodulators vs untreated
    * T and B cells vs untreated
* Figure 4 - PICRust bar plots
  * Plot treatment groups vs untreated (individually)
    * Immunodulators vs untreated
    * T and B cells vs untreated
* Supplemental figures
  * S1: Preliminary alpha and beta plots from PPMS vs SPMS
  * S2: Alpha and beta diversity from Aim 2 for individual drugs
  * Mention that Fingolimod shows lower diversity

## To-do
- [ ] Redo alpha and beta diversity for aim 2 with treatments grouped by mechanism
  * Do not need to run FDR for alpha/beta diversity for aim 2 (get rid of this)
  * PERMANOVA already gives adjusted p value
  * Not necessary for alpha because of the type of comparison
- [ ] Identify taxa from core microbiome
  * Earlier on in the code, a list is generated
Cross reference with OTU table and ASV IDs
- [ ] Annotate DESeq2 by only labeling the bacteria that are large effect
  * Make large effect >2.5 (do this anyways) and then label the large effect ones
  * Could talk about what these bacteria are in our discussion
  * ONLY Have new plots compare immunomodulators vs. untreated; t/b cells vs. untreated
  * Always have reference group be the control or untreated
- [ ] Change colours for aim2 (alpha and beta)
  * Change rituxan colour to #ffea00
- [ ] Run PICrust
  * Same comparisons as deseq2


