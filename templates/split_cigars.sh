#!/bin/bash
# https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads

echo "Loading GATK module..."
module load gatk

mkdir -p ${sample_id}
cd ${sample_id}
gatk SplitNCigarReads \
    -R "../${ref_path}" \
    -I "../${marked_dup_bam}" \
    -O ${sample_id}_split.bam

echo "COMPLETED!"
