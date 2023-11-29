#!/bin/bash

echo "\$(date)  Load modules..."
module load picard
module load gatk/3.8

echo "\$(date)  Create output folder..."
mkdir -p "${sample_id}"

echo "\$(date)  Create targets for realigning around indels..."
java -Xmx6g -jar \$gatk_dir/GenomeAnalysisTK.jar \
    -nt ${task.cpus} \
    -T RealignerTargetCreator \
    -R ${ref_path} \
    -I "${marked_dup_bam}" \
    -known "${indel_db1}" \
    -known "${indel_db2}" \
    -o "${sample_id}/${sample_id}_realigner.intervals"

touch "ok.txt"

echo "\$(date)  COMPLETED!"