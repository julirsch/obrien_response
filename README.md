# Confounding in _ex vivo_ models of Diamond Blackfan anemia

This repository contains all data and analyses performed to assess
potential confounding of the transcriptomic _ex vivo_ profiles
examined in [O'Brien _et al_](http://www.bloodjournal.org/content/early/2017/04/03/blood-2017-01-760462?sso-checked=true)
due to variability in cell type composition. Below is a roadmap
to our analysis framework as well as a synthesis of our results. 

## Raw Data

All raw microarray samples (.CEL files) can be found in the [data](data) subdirectory. 

- Study <-> GSE
- Study <-> GSE

## Code

- Including HTML

## Processed Data


## Overall

In brief, these analyses suggest that the _ex vivo_ samples produced by O'Brien _et al_
did not have identical cell composition, which we show can confound the association
reported in this study. Moreover, our results suggest that the effects uncovered in their
analyses can be attributable to cell type composition. The following figure summarizes our
findings from a re-analysis of their data. 

![Figure 1](media/Ulirsch_Figure1.png)
Figure 1. Evidence for confounding in microarray analysis of Diamond Blackfan anemia models.
(A) Diagram of hypothetical relationships between exposure, outcome, and putative confounding variable in the O’Brien et al. study. (B) Accurate deconvolution of early (CD235a-) and late (CD235a+) erythroid maturation stages by CIBERSORT in O’Brien et al. samples based upon normal erythroid maturation from GSE22552 (p = 0.000057 from likelihood ratio test). (C) Although similarly sorted for CD235a-, DBA samples are comprised of different mixtures of maturation stages than unaffected control samples. DBA due to RP or indeterminate (RP/I) samples are on average more mature (p = 0.017 from likelihood ratio test), whereas DBA due to GATA1 samples are less mature (p = 0.012 from likelihood ratio test). (D) Kernel density plots for heme biosynthesis, ribosome biogenesis, and curated (cur.) GATA1 target genes. Heme biosynthesis and ribosome biogenesis are significantly associated with erythroid maturation stage (p < 10^-10 for both by Kruskal-Wallis), whereas curated GATA1 target genes are only significant by more sensitive pairwise GSE analysis tests. Ribosome biogenesis was chosen as a more general measurement of “translational machinery” since snoRNAs were not measured in all microarrays. (E) Similar to (D), except for the three groups investigated in O’Brien et al. Differences between genotype (RP/I, GATA1, or unaffected) are much smaller than between stages (Kruskal-Wallis p > 0.5 for all comparisons, but significant by pairwise GSE anlaysis). (F) Similar GSE analysis to that reported by O’Brien et al. is shown. Synthetic stage-matched normal samples were created using the estimate from (B) of the percentage of each erythroid stage present. GSE analysis indicates that these synthetic normals have equivalent or stronger GSE results when compared directly to the original samples. Black bars indicate genes that are ranked according to expression differences between DBA (RP/I) and unaffected controls and normalized enrichment scores (NES) are reported. All data presented in this figure was RMA-normalized (see online analysis for SCAN-normalized). 
