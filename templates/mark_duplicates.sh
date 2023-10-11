#!/bin/bash


module load picard


java -jar $picard_dir/picard.jar MarkDuplicates \
    I=${run_id}_read_groups.bam \
    O=${run_id}_marked_dup.bam \
    M=${run_id}_marked_dup_metrics.txt
