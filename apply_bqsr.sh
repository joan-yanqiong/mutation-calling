#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name apply_bqsr_single	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load gatk

ref_path=$1
run_id=$2

gatk ApplyBQSR \
    -R ${ref_path} \
    -I ${run_id}_split.bam \
    --bqsr-recal-file ${run_id}_recal_data.table \
    -O ${run_id}_recal.bam