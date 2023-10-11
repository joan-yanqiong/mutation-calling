#!/bin/bash

# The command below was to retrieve the patient ID corresponding to the run ID from a csv file, but
# I think this is unnecessary. Any unique number for the read group will do, so the run_id should suffice.
# sample_id=$(cat ${csv_path} | grep ${run_id} | awk -F ',' '{print $2}')
module load picard/2.10.9

java -jar "${picard_dir}/picard.jar AddOrReplaceReadGroups" \
    I=${run_id}.bam \
    O=${run_id}_read_groups.bam \
    RGID=${i} \
    RGLB=lib_${i} \
    RGPL=ILLUMINA \
    RGPU=ILLUMINA_${run_id}_lib${i} \
    RGSM=${run_id}