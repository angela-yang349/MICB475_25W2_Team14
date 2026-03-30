# Chapter 6 - Aim 3

## Purpose:
Identify indicator taxa that are strongly associated with treatment status and specific treatment modalities in PMS patients

## Code:
Core Microbiome Analysis
 * [Aim 3a Code](Aim3_code.R) - original code
 * [Aim 3a Code 2 Groups](Aim3a_code2.R) - code modified to compare grouped treatments (GA/IFN/DMF vs ocrevus/fingolimod)
 * [Aim 3a Code 3 Groups](Aim3a_code_3_groups.R) - code modified to compare grouped treatments (DMF vs GA/IFN vs ocrevus/fingolimod)

Indicator Species Analysis
 * [Aim 3b Code](Aim3b_code.R) - original code containing ungrouped treatments and grouping by treatment_status (control vs treated vs untreated)
 * [Aim 3b Code 2 Groups](Aim3b_code_group2.R) - code modified to compare grouped treatments (GA/IFN/DMF vs ocrevus/fingolimod)
 * [Aim 3b Code 3 Groups](Aim3b_code_group3.R) - code modified to compare grouped treatments (DMF vs GA/IFN vs ocrevus/fingolimod)

## Methods:
* Core Microbiome Analysis
  * compare healthy controls vs untreated PMS vs treated PMS
  * compare healthy controls vs untreated PMS vs all the different PMS treatments
  * compare healthy controls vs untreated PMS vs immunomodulators (Glatiramer acetate, Interferon (IFN), Dimethyl fumarate) vs T/B cell treatments (Ocrevus (rituxan), Fingolimod)
  * compare healthy controls vs untreated PMS vs immunomodulators (Glatiramer acetate, Interferon (IFN)) vs T/B cell treatments (Ocrevus (rituxan), Fingolimod) vs DMF
    * abundance threshold of 0.001
    * prevalence threshold of 0.5
  * identified unique ASVs
* Indicator Species Analysis
  * compare healthy controls vs untreated PMS vs treated PMS
  * compare healthy controls vs untreated PMS vs all the different PMS treatments
  * compare healthy controls vs untreated PMS vs immunomodulators (Glatiramer acetate, Interferon (IFN), Dimethyl fumarate) vs T/B cell treatments (Ocrevus (rituxan), Fingolimod)
  * compare healthy controls vs untreated PMS vs immunomodulators (Glatiramer acetate, Interferon (IFN)) vs T/B cell treatments (Ocrevus (rituxan), Fingolimod) vs DMF

## Visualizations:
### Core microbiome analysis:

![MS Venn Diagram](ms_venn_general.png)

![MS Venn Diagram with Treatments](ms_venn_specific_treatments.png)

![MS Venn Diagram 2 Treatment Groups](venn_mechanism_1.png)

Lymphocyte:
|ASV |Taxa|
|----|----|
|5e22ac39af6385a9a5d7a57992e0d8ba |d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Rikenellaceae; g__Alistipes| 
|c073917b4d37d0b4f65cf4fb1b43619e|d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides|
|ce7260fe3f31c7fa47ae17dc33d2b8f7|d__Bacteria; p__Bacillota; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__Ruthenibacterium|
|73e5b567bf8428743df8b7f78c8ec1b0|d__Bacteria; p__Bacillota; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae; g__[Eubacterium]_siraeum_group|

Immunomodulators:
|ASV |Taxa|
|----|----|
|d52a0395431a1912fd52bda732cd2784|d__Bacteria; p__Actinomycetota; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Varibaculum|

Shared between lymphocyte and immunomodulators:
|ASV |Taxa|
|----|----|
|fc89b71ba14dc877e9e8527e6ed60391|d__Bacteria; p__Bacillota; c__Clostridia; o__Oscillospirales; f__Oscillospiraceae; g__Oscillibacter|

![MS Venn Diagram 3 Treatment Groups](venn_3_groups.png)

### Indicator species analysis:
Stat values must be 0.7 or above and indicator genuses should only belong to one group in order to be an indicator.

![ISA Ungrouped](ISA_results.csv)
Summary
* 1 indicator genus found with Glatiramer acetate: Peptococcus

![ISA Grouped By Treatment Status](ISA_results_treatment_status.csv)
Summary
* no indicators found

![ISA 2 Treatment Groups](ISA_results_group2.csv)
Summary
* no indicators found

![ISA 3 Treatment Groups](ISA_results_group3.csv)
Summary
* no indicators found

