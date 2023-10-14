#!/bin/bash

module load picard

mkdir -p "${sample_id}"

echo "1. Sorting bam file..."

java -jar \$picard_dir/picard.jar SortSam \
      I=${marked_dup_bam} \
      O=${sample_id}/${sample_id}_marked_dup_sorted_coord.bam \
      SORT_ORDER=coordinate

echo "2. Indexing bam file..."

module load samtools
samtools index ${sample_id}/${sample_id}_marked_dup_sorted_coord.bam

echo "Indexing done. Running GATK..."

module load gatk/3.8

echo "3. Create targets for realigning around indels..."
java -jar \$gatk_dir/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R ${ref_path} \
    -I "${sample_id}/${sample_id}_marked_dup_sorted_coord.bam" \
    -known $indel_db1 \
    -known $indel_db2 \
    -o "${sample_id}/${sample_id}_realigner.intervals"

echo "COMPLETED!"