#!/bin/bash

#SBATCH -n 41
#SBATCH -w ouga03
#SBATCH --mem 450000

module load i12g/R/3.4.2-Bioc3.5

Rscript src/r/processing/get_counts_batch_3.R -n 40