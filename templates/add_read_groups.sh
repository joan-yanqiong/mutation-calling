#!/bin/bash
echo "\$(date)  Loading picard module..."
module load picard/2.10.9

echo "\$(date)  Creating directory for output..."
mkdir -p "${sample_id}"

echo "\$(date)  Check sample type..."
if [[ ${sample_type} == "tumor" ]]
then
    lib_id=${ix}
elif [[ ${sample_type} == "normal" ]]
then
    lib_id=1
fi

echo "\$(date)  Adding read groups to ${sample_id}..."
java -jar \${picard_dir}/picard.jar AddOrReplaceReadGroups \
    I="${mapped_bam}" \
    O="${sample_id}/${sample_id}_RG.bam" \
    RGID=\${lib_id} \
    RGLB=lib_\${lib_id} \
    RGPL=ILLUMINA \
    RGPU=ILLUMINA_${sample_id}_lib\${lib_id} \
    RGSM=${sample_id}

touch "ok.txt"
echo "\$(date)  COMPLETED!"