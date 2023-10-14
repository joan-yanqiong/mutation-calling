#!/bin/bash
echo "Loading picard module..."
module load picard/2.10.9

echo "Creating directory for output..."
mkdir -p "${sample_id}"

echo "Adding read groups to ${sample_id}..."
java -jar \${picard_dir}/picard.jar AddOrReplaceReadGroups \
    I="${mapped_bam}" \
    O="${sample_id}/${sample_id}_Aligned.sortedByCoord.out_read_groups.bam" \
    RGID=${ix} \
    RGLB=lib_${ix} \
    RGPL=ILLUMINA \
    RGPU=ILLUMINA_${sample_id}_lib${ix} \
    RGSM=${sample_id}

echo "COMPLETED!"