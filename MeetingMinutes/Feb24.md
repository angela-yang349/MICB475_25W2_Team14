# February 24, 2026
## Agenda
  * Completed work:
    * completed proposal
  * Upcoming work:
    * PPMS and SPMS beta diversity
    * Aim 1: alpha diversity of PMS treated, PMS untreated, and healthy controls
  * Questions:
    * how should we be using the healthy controls?
    * can we go over the statistical methods for each step?
    * is there a way to correct for confounding variables like bmi, diet, depression, etc?
    * should we include a heatmap for ISA visualization in addition to the table?

## Discussion
### Ritu’s updates & suggestions
 * Don’t wait to complete the earlier analysis before starting. Once you have the phyloseq object, you can begin the later analysis 
  * Run parallel analyses 
* March 13/23 evelyn will join for 2 meetings so get most of the analysis done by that time – will be easier to discuss & can go over results/have a story 
 * Avoid waiting another week to get results 
* Start Google Slides; add plots onto that so all plots will be in the same place and easy to share with Ritu & Evelyn
 * Make them prettier later on for the presentation
* Next week: getting proposals back
 * Will get to revise and incorporate the revisions to get some points back
 * Can go over revisions in next weeks meeting (march 3rd (?))

### Currently
* Will be running diversity analyses 
 * Use healthy controls as a baseline to see whether PMS treated groups and healthy microbiome groups converge? 
 * 3-way analysis – all 3 together 
  * Analysis is the same
  * Just have another group compared to the treated and untreated group
  * Is the treated group becoming more similar to the healthy one?
  * DESEQ is a pairwise comparison
   * Treated vs untreated
   * Treated vs healthy
  * Venn diagram from core microbiome
   *For now, try to include the healthy group (include everything)
* Statistical methods
  * Indicator species analysis
   * May need to adjust abundance and and prevalence
  * Alpha diversity tells you if there’s significance in an individual group. Shannon diversity – within samples within one category
   * Kruskal-Wallis statistical test helps to compare the three groups (treated vs untreated vs healthy) tells you if a group is more diverse than the other?
   * Alpha first
   * Do shannon + kruskal-wallis for all alpha diversity analyses
  * Beta diversity is a direct comparison between the three 
   * Beta second
   * First time treated vs untreated, second time individual treatments 
   * Can do 2 sets of alpha diversity
   * Then proceed with ones that show a difference
   * Maybe binn certain treatments if they’re similar enough so fewer conditions and more samples within conditions
* Removing confounding variables
  * If you see an extreme sample (outlier), remove that one sample manually 
  * Eg., look at common known factors that can change the microbiome (age), skim through to see if any outliers in that parameter compared to others. If there’s something extreme, remove that
  * Keep samples that have overlaps in those groups to remove the confounding effect of eg. age
  * Keep it as it as unless there’s an outlier
  * Mention regions/treatment types in the discussion
* Heatmap for ISA visualization + table
  * Table needed for sure (gives # of indicator species for each category)
  * If lots of indicator species, make sure to mention the # for the condition
  * If not enough time, stick with the table, otherwise can have heatmap 

### Next steps:
Isabella: look at beta between SSMS and PPMS and see if we bin into PMS
Angela & Sadaf: Alpha and beta diversity of treated and untreated
