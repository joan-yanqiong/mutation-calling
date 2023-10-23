#!/bin/bash

# Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

echo "\$(date)  Loading picard module..."
module load picard

echo "\$(date)  Creating directory for output..."
mkdir -p ${sample_id}

echo "\$(date)  Marking duplicates..."
java -jar \${picard_dir}/picard.jar MarkDuplicates \
    I="${read_groups_file}" \
    O="${sample_id}/${sample_id}_marked_dup.bam" \
    M="${sample_id}/${sample_id}_marked_dup_metrics.txt"

echo "\$(date)  COMPLETED!"