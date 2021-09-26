#!/bin/bash

# run like
# sh 00_run_slim_sim.sh slim_sim/recomb_pop_evol_one_rate.slim data/slim_output/output_onerate_vcf.vcf test_onerate

slim_sim=$1
seq_len=10000000
slim_vcf=$2
out_fasta_slug=$3

slim -l $1

Rscript 01_slim_vcf_to_fasta.R $seq_len $slim_vcf $out_fasta_slug
