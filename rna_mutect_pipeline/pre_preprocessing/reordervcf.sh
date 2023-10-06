#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name reorder_vcf		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"

module load picard

java -jar $picard_dir/picard.jar SortVcf \
      I="${prefix}/Reference/Cosmic_GenomeScreensMutant_Normal_v98_GRCh37.vcf" \
      O="${prefix}/Reference/Cosmic_GenomeScreensMutant_Normal_v98_GRCh37_sorted.vcf" \
      SD="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.dict"