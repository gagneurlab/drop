---
title: Detection of RNA outliers Pipeline
---

## Mono-Allelic Expression Module 

Go to the MAE tab and click on `Datasets.R`. Then, select the desired 
analysis group. 
It contains a cascade plot with the number of allelic counts that passed the filter 
and were mono-allelically expressed, and the results table.

Go to the QC tab and click on `Datasets.R`. Then, select the desired 
analysis group. 
It contains a histogram with the percentage of DNA - RNA variants, and tables containing
the mismatches, if any.

MAE is computed using a negative binomial test described in [Kremer et al, Nat Commun 2017](https://www.nature.com/articles/ncomms15824).
