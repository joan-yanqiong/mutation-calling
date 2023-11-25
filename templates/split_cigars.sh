#!/bin/bash
# https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads

echo "\$(date)  Loading GATK module..."
module load gatk

echo "\$(date)  Creating directory for output..."
mkdir -p ${sample_id}

echo "\$(date)  Entering created directory..."
cd ${sample_id}

echo "\$(date)  Splitting reads with N in CIGAR string..."
gatk SplitNCigarReads --java-options -Xmx4g \
    -R "../${ref_path}" \
    -I "../${marked_dup_bam}" \
    -O ${sample_id}_split.bam

touch "ok.txt"

echo "\$(date)  COMPLETED!"
