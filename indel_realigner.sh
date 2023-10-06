#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-20:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name indel_realigner	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

run_id=$1
ref_path=$2
indel_db1=$3
indel_db2=$4

module load gatk/3.8

java -jar $gatk_dir/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ${ref_path} \
    -I ${run_id}_marked_dup_sorted_coord.bam \
    -known $indel_db1 \
    -known $indel_db2 \
    -targetIntervals ${run_id}_realigner.intervals \
    -o ${run_id}_realigned.bam