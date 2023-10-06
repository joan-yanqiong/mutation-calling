#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name bqsr_table_single	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load gatk

ref_path=$1
dbSNP_vcf=$2
run_id=$3

gatk BaseRecalibrator \
    -R ${ref_path} \
    -I ${filename}_split.bam \
    --known-sites ${dbSNP_vcf} \
    -O ${run_id}_recal_data.table