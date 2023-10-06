#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name bqsr_normal_table_single	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load gatk

run_id=$1
ref_path=$2
dbSNP_vcf=$3

gatk BaseRecalibrator \
    -R ${ref_path} \
    -I ${run_id}_realigned.bam \
    --known-sites ${dbSNP_vcf} \
    -O ${run_id}_recal_data.table