#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=40G 
#SBATCH --job-name mark_duplicates_and_split_cigars		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load picard

run_id=$1

java -jar $picard_dir/picard.jar MarkDuplicates \
    I=${run_id}_read_groups.bam \
    O=${run_id}_marked_dup.bam \
    M=${run_id}_marked_dup_metrics.txt
