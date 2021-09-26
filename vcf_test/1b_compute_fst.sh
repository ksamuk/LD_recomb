#!/bin/bash

#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH -p noor
#SBATCH -J fst-%j
#SBATCH --mail-type=END,FAIL
#SBATCH -o tmp/fst-%j.out

# compute pairwise fst between populations over course of slim simulations

# get list of param combo folders

ls data/slim > tmp/parm_combo_folders.txt

while read paramfolder; do
	ls $paramfolder > tmp/repfolder.txt
	while read repfolder; do
		ls $simfolder/$repfolder/*.gz.gz | sed 's/.*\/p[12]_//g' | sed 's/_.*//g' | uniq > tmp/fst_genlist.txt
		while read gen; do
			vcf1="$simfolder/$repfolder/p1_${gen}_tmp.vcf.gz.gz"
			vcf2="$simfolder/$repfolder/p2_${gen}_tmp.vcf.gz.gz"
			bcftools merge $vcf1 $vcf2 --force-samples > tmp/fst_tmp.vcf
			fst=$(vcftools --vcf fst_tmp.vcf --fst-window-size 100000 --weir-fst-pop meta/fst_pops_1.txt --weir-fst-pop meta/fst_pops_2.txt -c | awk 'FNR == 2 {print $5}')
			echo "$paramfolder $repfolder $vcf1 $vcf2 $fst" >> data/fst_comparisons.txt
		done < tmp/fst_genlist.txt
	done < tmp/repfolder.txt
done < tmp/parm_combo_folders.txt