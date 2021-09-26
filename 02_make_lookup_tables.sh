#!/bin/bash

#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH -p noor
#SBATCH -J make_tab-%j
#SBATCH --mail-type=END,FAIL
#SBATCH -o tmp/make_tab-%j.out

# make a lookup table for pyrho

mkdir -p lookup_tables

echo "Creating lookup table..."

# create separate lookup tables for each generation post divergence
for gen in `seq 1721 1721 51630`
do 

echo "CREATING PYRHO LOOKUP TABLE FOR GENERATION $gen"

pyrho make_table --samplesize 50 --numthreads 4 \
--approx  --moran_pop_size 75 \
--mu 5.49e-09 \
--popsizes 860300,1720600 \
--epochtimes $gen \
--outfile lookup_tables/slim_sim_lookup_splitAware_gen_${gen}.hdf

done

extragen="172100 1721000"
for gen in $extragen
do 

echo "CREATING PYRHO LOOKUP TABLE FOR GENERATION $gen"

pyrho make_table --samplesize 50 --numthreads 4 \
--approx  --moran_pop_size 75 \
--mu 5.49e-09 \
--popsizes 860300,1720600 \
--epochtimes $gen \
--outfile lookup_tables/slim_sim_lookup_splitAware_gen_${gen}.hdf

done

