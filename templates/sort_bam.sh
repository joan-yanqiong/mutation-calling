#!/bin/bash

echo "\$(date)    Loading picard module..."
module load picard

echo "\$(date)    Creating directory for output..."
mkdir -p ${sample_id}

echo "\$(date)    Sorting ${sample_id}..."
java -Xmx3g -jar \$picard_dir/picard.jar SortSam \
      I=${read_groups_sam} \
      O=${sample_id}/${sample_id}_sortedBy_${sort_order}.bam \
      SORT_ORDER=${sort_order}
touch "ok.txt"
echo "\$(date)    COMPLETED!"