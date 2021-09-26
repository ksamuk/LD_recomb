#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=20
#SBATCH -p noor
#SBATCH -J slim_RE
#SBATCH -o tmp/slim_RE-%j.out


Rscript 01_run_recomb_evolve_realistic_rates_gene_flow.R

