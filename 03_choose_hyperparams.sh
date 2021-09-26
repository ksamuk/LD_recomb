#!/bin/bash

# make a lookup table for pyrho

mkdir -p hyperparams

pyrho hyperparam --samplesize 50 \
--tablefile lookup_tables/slim_sim_lookup.hdf \
--mu 1.4e-09 \
--blockpenalty 50,100 \
--windowsize 25,50 \
--num_sims 3 \
--ploidy 2 \
--popsizes 1250000,1250000 \
--epochtimes 10000 \
--outfile hyperparams/slim_sim_hyperparams.txt