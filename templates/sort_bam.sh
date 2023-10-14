#!/bin/bash

echo "Loading picard module..."
module load picard

echo "Creating directory for output..."
mkdir -p ${sample_id}

echo "Sorting ${sample_id}..."
java -jar \$picard_dir/picard.jar SortSam \
      I=${read_groups_bam} \
      O=${sample_id}/${sample_id}_sorted.bam \
      SORT_ORDER=queryname

echo "COMPLETED!"