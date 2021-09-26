#!/bin/bash
#SBATCH -J test_array
#SBATCH -p noor
#SBATCH -t 0-2:00
#SBATCH --mem 8000
#SBATCH -o tmp/array_test-%A_%a.out.out

# launches a job array of recombination rate evolution simulations
# used by Rscript launcher to launch arrays for each paramter combination
# individual arrays can be run like:
# sbatch --array=0-1 scripts/run_slim_sim_array.sh 1000 1000 0.00000175 0.00003075 1 1 100000 1000 0.5 data/slim_output/1000_1000_0.00000175_0.00003075_1_1_100000_1000_0.5

echo "RUNNING SLIM SIMULATION WITH ARRAY INDEX ${SLURM_ARRAY_TASK_ID}"
echo "input was : ${1}"
