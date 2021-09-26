#!/bin/bash
#SBATCH -J slim_array-%j
#SBATCH -t 0-2:00
#SBATCH --mem 8000
#SBATCH -o tmp/slim_RE/slim_RE-%A_%a.out

# launches a job array of recombination rate evolution simulations
# used by Rscript launcher to launch arrays for each paramter combination
# individual arrays can be run like:
# sbatch --array=0-3 scripts/run_slim_sim_array.sh p1si=1000 p2si=1000 mu=0.00000175 rr=0.00003075 rrc=1.0 rrcd=1 chrlen=100000 samp=1000 mig=0.5 outfold=data/slim_output/1000_1000_0.00000175_0.00003075_1_1_100000_1000_0.5/

echo "RUNNING SLIM SIMULATION WITH ARRAY INDEX ${SLURM_ARRAY_TASK_ID}"

# parse named arguments
# the array syntax does something weird to numbered arguments
# so this is the workaround
for ARGUMENT in "$@"
do

    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=)   

    case "$KEY" in
            p1si)    p1si=${VALUE} ;;
            p2prop)    p2prop=${VALUE} ;;
            rrc)     rrc=${VALUE} ;;
            rrcd)    rrcd=${VALUE} ;;
            varrr)    varrr=${VALUE} ;;
            chrlen)  chrlen=${VALUE} ;;
            samp)    samp=${VALUE} ;;
            mig12)     mig12=${VALUE} ;;
            mig21)     mig21=${VALUE} ;;
            outfold) outfold=${VALUE} ;;
            split)   split=${VALUE} ;;
            gfgen)   gfgen=${VALUE} ;;
            *)   
    esac    

done

# make the base  output directory 
mkdir -p $outfold

slim -d p1_size=$p1si \
-d p2_prop=$p2prop \
-d recomb_rate_cha=$rrc \
-d recomb_rate_change_dir=$rrcd \
-d variable_recomb=$varrr \
-d chr_length=$chrlen \
-d samp_interval=$samp \
-d mig_rate_m12=$mig12 \
-d mig_rate_m21=$mig21 \
-d split_generation=$split \
-d gene_flow_generation=$gfgen \
-d out_folder_ba=\'${outfold}\' \
slim_sim/recomb_rate_evolve.slim