#!/bin/bash

#SBATCH --mem=220000
#SBATCH --cpus-per-task=22
#SBATCH -p noor
#SBATCH -J to_fasta
#SBATCH -o slurm-to_fasta.out

Rscript 02_post_process_slim_to_fasta.R

