#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name mark_duplicates_trial		
#SBATCH -o %x-%j.out

module load picard

run_id=$1

echo "Running MarkDuplicates on ${run_id}."

java -jar $picard_dir/picard.jar MarkDuplicates \
    I=${run_id}_sorted.bam \
    O=${run_id}_marked_dup.bam \
    M=${run_id}_marked_dup_metrics.txt