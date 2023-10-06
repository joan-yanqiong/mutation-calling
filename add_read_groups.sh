#!/bin/bash

#SBATCH -p all
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=4G 
#SBATCH --job-name read_groups		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

i=$1
run_id=$2

# The command below was to retrieve the patient ID corresponding to the run ID from a csv file, but
# I think this is unnecessary. Any unique number for the read group will do, so the run_id should suffice.
# sample_id=$(cat ${csv_path} | grep ${run_id} | awk -F ',' '{print $2}')

java -jar $picard_dir/picard.jar AddOrReplaceReadGroups \
    I=${run_id}.bam \
    O=${run_id}_read_groups.bam \
    RGID=${i} \
    RGLB=lib_${i} \
    RGPL=ILLUMINA \
    RGPU=ILLUMINA_${run_id}_lib${i} \
    RGSM=${run_id}