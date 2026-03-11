# March 10th, 2026
## Agenda
* Completed work:
  * Removal of RRMS household controls
  * Core microbiome
  * Indicator species analyses
  
* Upcoming work:
  * DEseq2
  * PICRUSt2

* Questions:
  * do we need to make a bar plot for our ISA results?
  * how do we upload the ISA table to GitHub?

## Conclusions
*  Do research (find papers) that can justify our treatment groupings by mechanism of action
*  Run core microbiome and ISA again with new groupings and any relevant changes

## Discussion
* Core microbiome:
  * Run code again but be less stringent (But we are already being very lenient)
  * Make objects for grouping treatments by mechanism
* Indicator Species Analysis:
  * Different numbers of significant taxa between runs is expected because multipatt() uses permutations.
    * Use set.seed() to make it reproducible
  * seed was set to 100 arbitrarily. Got 7 genuses that are indicative
  * just having the ISA table is fine for visualizing the results
  * Run ISA again with treatments grouped all together
  * Run ISA once more with treatments grouped by mechanism
* Types of treatments for grouping:
  * source: https://www1.racgp.org.au/ajgp/2022/april/multiple-sclerosis-diagnosis-therapy-and-prognosis
  * Dimethyl fumarate (6 samples) ORAL → Antioxidative, immunomodulatory and anti-inflammatory effects via activation of the nuclear factor (erythroid-derived 2)-like 2 transcriptional pathway
  * Fingolimod (3 samples) → ORAL Sphingosine 1-phosphate receptor modulators; prevent lymphocyte trafficking through the lymph node and cause a reversible lymphopaenia
  * Glatiramer acetate (2 samples) → SUBCUTANEOUS Synthetic polypeptide; possibly blocks presentation of myelin antigens to T lymphocytes
  * Interferon (8 samples) → SUBCUTANEOUS Immunoregulatory actions including antagonism of interferon gamma, reduction of cytokine release and augmentation of suppressor T cell function
  * Anti-CD20 monoclonal antibodies = ocrevus (rituxan) (12 samples) → INFUSED Humanised monoclonal antibody to CD20, depletes B lymphocytes
  * Natalizumab (0) → Recombinant humanised monoclonal antibody to α4-integrins, inhibiting leucocyte migration from blood to CNS (RRMS only)
    * No longer considered in our dataset as it was only used for RRMS treatment
* Groupings:
  * Mechanism of action:
    * Immunomodulators (moderate efficacy) = Glatiramer acetate, Interferon (IFN), Dimethyl fumarate
    * T cell/B cell (high efficacy) = Ocrevus (rituxan), Fingolimod
  * Mechanism of action (stringent):
    * Antioxidant =  Glatiramer acetate
    * Immunomodulators = Interferon (IFN), Dimethyl fumarate
    * T cell/B cell = Ocrevus (rituxan), Fingolimod
  * Treatment delivery:
    * Oral: Dimethyl fumarate, Fingolimod
    * Subcutaneous: Glatiramer acetate, Interferon (IFN)
    * Infused: Ocrevus (rituxan)

## To-do
- [ ] Make sure to set.seed (100) for ISA
- [ ] Run an ISA with all treatments in one category and then grouped by mechanism
- [ ] Run core microbiome again with treatments grouped by mechanism
- [ ] Try to remove zero labels in the Venn diagrams for clarity
- [ ] Save ms_rare with RRMS controls removed into a new RData file
- [ ] Do research (find papers) that can justify our treatment groupings by mechanism of action

