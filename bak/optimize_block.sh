#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=20
#SBATCH -p noor

Rscript 02_optimize_block_window_ldhelmet.R $1
