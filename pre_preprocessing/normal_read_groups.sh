#!/bin/bash

#SBATCH -p all
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=4G 
#SBATCH --job-name normal_read_groups		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

run_id=$1
i=$2

java -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
    I=${run_id}_mapped.bam \
    O=${run_id}_read_groups.bam \
    RGID=${i} \
    RGLB=lib_${i} \
    RGPL=ILLUMINA \
    RGPU=ILLUMINA_${run_id}_lib${i} \
    RGSM=${run_id}