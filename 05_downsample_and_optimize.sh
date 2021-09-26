#!/bin/bash

# downsample and run the optimization pupline

# downsample vcfs to match lookup tables
vcf=$1
vcfout=$(echo $vcf | sed 's/\.vcf/_n25.vcf/g')

gunzip -c $vcf | bgzip > $vcfout

# determine generation for the lookup table
gen=$(echo $vcf | sed 's/.*_//g' | sed 's/\.vcf.gz//g')

# lookup table generations are relative to divergence (t=20000)
let gen=$(gen-20000)

# optimize with pyrho

mkdir -p optimize/split_aware/$2/$3
mkdir -p optimize/split_unaware/$2/$3

out=$(echo $vcfout | sed 's/.*\///g' | sed 's/.vcf.gz/_recomb.rmap/g')

pyrho optimize --vcffile $vcfout \
--windowsize 100 \
--blockpenalty 1000 \
--tablefile lookup_tables/slim_sim_lookup_splitAware_gen_${gen}.hdf \
--ploidy 2 \
--numthreads 4 \
--outfile optimize/split_aware/$2/$3/$out

rm $vcfout
