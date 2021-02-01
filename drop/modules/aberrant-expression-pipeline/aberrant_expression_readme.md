---
title: Detection of RNA outliers Pipeline
---

## Aberrant Expression Module 

Go to the Counting tab and click on `Datasets.R`. Then, select the desired 
analysis group. 
It contains different plots like the number of reads counted, size factors, 
expression of each gene in FPKM and the number of expressed genes. 

Go to the Outrider tab and click on `Datasets.R`. Then, select the desired 
analysis group. 
It contains different plots like the steps towards finding the optimal encoding 
dimension, heatmaps of the sample covariation before and after correction, 
biological coefficient of variation; and tables containing aberrant samples and 
results table.

Expression outliers are computed using [OUTRIDER](https://www.cell.com/ajhg/fulltext/S0002-9297(18)30401-4) 
(Brechtmann, Mertes et al., AJHG, 2018).
