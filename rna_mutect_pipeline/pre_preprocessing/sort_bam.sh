#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=20G
#SBATCH --job-name sort_bam
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

run_id=$1

module load picard

java -jar $picard_dir/picard.jar SortSam \
      I=${run_id}_read_groups.bam \
      O=${run_id}_sorted.bam \
      SORT_ORDER=queryname
