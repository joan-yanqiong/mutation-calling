#!/bin/bash

echo "Loading picard module..."
module load picard

echo "Creating directory for output..."
mkdir -p ${sample_id}

echo "Sorting ${sample_id}..."
java -jar \$picard_dir/picard.jar SortSam \
      I=${read_groups_sam} \
      O=${sample_id}/${sample_id}_${suffix}.bam \
      SORT_ORDER=${sort_order}

echo "COMPLETED!"