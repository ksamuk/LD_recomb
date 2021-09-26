#!/bin/bash

#SBATCH -p noor
#SBATCH --cpus-per-task=8
#SBATCH -J process_pyrho
#SBATCH --mail-user=ksamuk@gmail.com
#SBATCH --mail-type=END,FAIL
#SBATCH -o tmp/process_pyrho-%j.out

Rscript 09_process_pyrho_files.R

