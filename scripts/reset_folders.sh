#!/bin/bash


rm -r empty dir
mkdir -p empty_dir
rsync -a --delete empty_dir/    data/slim_output

rm -r empty_dir
rm -r data/slim_output
rm -r tmp/slim_RE

rm tmp/*

mkdir -p data/slim_output
mkdir -p tmp/slim_RE

sh 00_run_sims_realistic_gene_flow.sh
