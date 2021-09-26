#!/bin/bash

# make a lookup table for pyrho

mkdir -p optimize/$1/$2

vcf=$3
out=$(echo $vcf | sed 's/.*\///g' | sed 's/.vcf.gz/_recomb.rmap/g')

pyrho optimize --vcffile $vcf \
--windowsize 100 \
--blockpenalty 1000 \
--tablefile lookup_tables/slim_sim_lookup_n50_splitAware.hdf \
--ploidy 2 \
--numthreads 1 \
--outfile optimize/$1/$2/$out

rm $vcf