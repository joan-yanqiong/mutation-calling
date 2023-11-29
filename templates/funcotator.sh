#!/bin/sh

echo "\$(date)   Loading modules..."
module load gatk

echo "\$(date)  Creating output directory..."
mkdir -p ${sample_id}


# output=\$(${params.gatk_path} CountVariants -V ${input_vcf})
# n_variants=\$(echo \$output | grep -oP '(?<=Tool returned: )[0-9]+')

# if [[ \$n_variants -eq 0 ]] ; then
#     echo "No variants found!"
#     touch "ok.txt"
#     exit 0
# else
echo "\$(date)  Use Funcotator for annotation..."
${params.gatk_path} Funcotator \
    --variant ${input_vcf} \
    --reference ${ref_path} \
    --ref-version hg19 \
    --data-sources-path ${params.func_data_path} \
    --output ${sample_id}/${sample_id}_${suffix}.funcotated.maf \
    --output-file-format MAF
# fi

touch "ok.txt"

echo "\$(date)   COMPLETED!"