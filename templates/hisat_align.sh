#!/bin/bash

echo "\$(date)\tLoading hisat2 module..."
module load hisat2/2.0.4

echo "\$(date)\tSample type: ${sample_type}"
if [[ ${sample_type} == "normal" ]]
then
    echo "\$(date)\tPair ID: ${pair_id}"
    output_folder=${normal_id}
    output_file_prefix=${pair_id}
elif [[ ${sample_type} == "tumor" ]]
then
    echo "\$(date)\tTumor ID: ${tumor_id}"
    output_folder=${tumor_id}
    output_file_prefix=${tumor_id}
fi

echo "\$(date)\tCreate output folder..."
mkdir -p \${output_folder}

echo "\$(date)\tAligning reads..."
# Parameters from Supp Table 12
hisat2 -p "${task.cpus}" -k 5 \
    --min-intronlen 20 \
    --max-intronlen 500000 \
    -X 800 \
    -I 0 \
    -x "${index_dir}/${index_dir}" \
    -1 "${fastq1}" \
    -2 "${fastq2}" \
    -S "\${output_folder}/\${output_file_prefix}_aligned_hisat2.sam" \
    --rg-id=1 \
    --rg "SM:\${output_file_prefix}\\tPU:ILLUMINA_\${output_file_prefix}_lib1" \
    --rg "LB:lib_1" \
    --rg "PL:ILLUMINA"
echo "\$(date)\tCompleted"