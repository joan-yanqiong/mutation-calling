#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name mark_duplicates_and_split_cigars		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

ref_path=$1
run_id=$2

gatk SplitNCigarReads \
    -R ${ref_path} \
    -I ${run_id}_marked_dup.bam \
    -O ${run_id}_split.bam