# March 3rd, 2026
## Agenda
* Completed work:
  * compared PPMS and SPMS samples through beta-diversity (aim 1)
  * PPMS and SPMS samples were binned together into PMS patients for all subsequent analyses (aim 1)
  * compared treated PMS vs. untreated PMS vs. healthy controls through alpha diversity (aim 1)
  * compared treated PMS vs. untreated PMS vs. healthy controls through beta diversity (aim 1)
  * compared each type of treatment PMS vs. untreated PMS vs. healthy controls through alpha diversity (aim 2)
  * compared each type of treatment PMS vs. untreated PMS vs. healthy controls through beta diversity (aim 2)
* Upcoming work:
  * core microbiome analysis (aim 3)
  * indicator species analysis (aim 3)
* Questions:
  * project 2 proposal revisions?
  * conducted pairwise PERMANOVA comparisons for beta diversity (aim 1 & aim 2) - did we do them correctly?
  * differences between bonferroni and FDR correction? is this necessary for PERMANOVA?
  * is pairwise Wilcoxon tests with FDR correction needed for alpha diversity results?
  * are individual ANOVA tests for alpha diversity only run is significant is detected using the Kruskal-Wallis test?

## Conclusions
* Remove RRMS household controls from the data and assign it to a new object
* Re-run the diversity analyses again using the data with the removed RRMS controls
* Do proposal corrections

## Discussion
* FDR correction is a good call
* Adjust abundance and prevalence thresholds for core microbiome to get a few results
  * Just for representation purposes
* Make a table with the number of indicator species found for each condition if there are many
* For the N/A in indicator species, go back to rep_seqs and BLAST the ASVs to get a probable species
* Remove RRMS household controls
* Maybe pair treatment conditions since they are low in number

## To-do
- [ ] proposal corrections
- [ ] Re-run the diversity analyses again using the data with the removed RRMS controls
