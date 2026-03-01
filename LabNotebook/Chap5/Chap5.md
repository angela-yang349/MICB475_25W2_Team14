# Chapter 5 - Aim 2

## Purpose:
To evaluate whether distinct disease-modifying therapies are significantly associated with differences in gut microbiome profiles among PMS patients.

## Code:
[Aim 2 code](/Aim2_code.R)

## Visualizations:

### Alpha diversity untreated PMS vs. each treatment type (with controls)

![treatment_types_shannon_plot](treatment_types_shannon_plot.png)

No significant differences in alpha diversity (Shannon Index) were detected across individual treatment types, including untreated and healthy control participants (Kruskal–Wallis, χ² = 8.89, df = 6, p = 0.180), indicating that within-sample microbial diversity did not vary significantly by treatment type.


### Beta diversity untreated PMS vs. each treatment type (with controls)

![treatment_types_wunifrac_pcoa_with_controls](treatment_types_wunifrac_pcoa_with_controls.png)


### Beta diversity untreated PMS vs. each treatment type (no controls)

![treatment_types_wunifrac_pcoa_pms_only](treatment_types_wunifrac_pcoa_pms_only.png)

Pairwise PERMANOVA tests revealed that none of the treatment comparisons showed statistically significant microbiome differences after FDR correction for multiple testing (all p_FDR > 0.05), despite several significant results before correction (fingolimod vs healthy controls: p = 0.026; glatiramer acetate vs healthy controls: p = 0.033).
Effect sizes were consistently small across all comparisons (R² < 0.03), with treatment type explaining less than 3% of microbiome variance in any pairwise comparison.

Pairwise PERMANOVA results: [treatment_types_pairwise_permanova](/treatment_types_pairwise_permanova.csv)
