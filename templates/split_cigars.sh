#!/bin/bash

gatk SplitNCigarReads \
    -R ${ref_path} \
    -I ${run_id}_marked_dup.bam \
    -O ${run_id}_split.bam