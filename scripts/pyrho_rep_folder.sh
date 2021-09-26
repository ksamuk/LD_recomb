#!/bin/bash

#SBATCH --mem=10000
#SBATCH --cpus-per-task=4
#SBATCH -p noor
#SBATCH -J optimfolder-%j
#SBATCH -o tmp/optim/optimfolder-%j.out

rootfolder=$1
simfolder=$2

#echo $repfolder

# list all the folder in the simulation folder
ls $rootfolder/$simfolder > tmp/rep_folder_list_${simfolder}.txt

# build a list of generations from the simulation to target 
# (19998 & 19999 are "definitely before split & gene flow")
genlist=(`counter=50000; seq -s\| $counter -1000 20000`)
genlist='19998|19999|'$genlist

while read repfolder
do

    ls $rootfolder/$simfolder/$repfolder/*.vcf.gz | grep -E $genlist > tmp/vcf_list_${simfolder}.txt

    while read vcf
    do
       	
         vcfout=$(echo $vcf | sed 's/\.vcf/_n25.vcf/g')
         optimout=$(echo $vcfout | sed 's/.*\///g' | sed 's/.vcf.gz/_recomb.rmap/g')
         			
         echo "Running pyrho on $vcfout..."
         			
         if [ -f "optimize/$simfolder/$repfolder/$optimout" ]; then
           echo "optimize/$simfolder/$repfolder/$optimout already exists, skipping..."
         fi
         
         if [ ! -f "optimize/$simfolder/$repfolder/$optimout" ]; then

           sh 05_downsample_and_optimize.sh $vcf $simfolder $repfolder
           
           #echo "05_downsample_and_optimize.sh $vcf $simfolder $repfolder"
           rm -f $vcfout
         fi
    
    done < tmp/vcf_list_${simfolder}.txt
  
done < tmp/rep_folder_list_${simfolder}.txt