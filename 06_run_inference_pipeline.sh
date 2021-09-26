#!/bin/bash

# run the pyrho pipeline for multiple vcfs


mkdir -p optimize
mkdir -p tmp/optim
rootfolder="/work/kms173/recomb_pop_diff/data/slim_output"
#simfolder="1000_0.5_0.0_-1_0_100000_1000_0.001_0.001_20000_30000"

ls $rootfolder  > tmp/sim_folder_list.txt

while read simfolder
do
	
		sbatch scripts/pyrho_rep_folder.sh  $rootfolder $simfolder
   
done < tmp/sim_folder_list.txt
