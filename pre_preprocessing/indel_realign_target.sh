#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-20:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name index_and_indel_realign_target	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

run_id=$1
ref_path=$2
indel_db1=$3
indel_db2=$4

# module load picard

# echo "Sorting bam file..."

# java -jar $picard_dir/picard.jar SortSam \
#       I=${run_id}_marked_dup.bam \
#       O=${run_id}_marked_dup_sorted_coord.bam \
#       SORT_ORDER=coordinate

echo "Indexing bam file..."

module load samtools
samtools index ${run_id}_marked_dup_sorted_coord.bam

echo "Indexing done. Running GATK..."

module load gatk/3.8

java -jar $gatk_dir/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ${ref_path} \
    -I ${run_id}_marked_dup_sorted_coord.bam \
    -known $indel_db1 \
    -known $indel_db2 \
    -o ${run_id}_realigner.intervals