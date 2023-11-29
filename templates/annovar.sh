#!/usr/bin/env bash
echo "\$(date)  Load module..."
module load annovar
module load gatk
echo "\$(date)   Create output directory..."
mkdir -p ${sample_id}


# output=\$(${params.gatk_path} CountVariants -V ${mutect_vcf})
# n_variants=\$(echo \$output | grep -oP '(?<=Tool returned: )[0-9]+')

# if [[ \$n_variants -eq 0 ]] ; then
#     echo "No variants found!"
#     touch "ok.txt"
#     exit 0
# else
    echo "\$(date)   Annotating variants..."
    perl ${params.annovar}/table_annovar.pl \
                    -vcfinput "${mutect_vcf}" \
                    ${human_db} \
                    -buildver hg19 \
                    -remove \
                    -protocol refGene \
                    -operation g \
                    -nastring . \
                    -polish \
                    --otherinfo \
                    -out "${sample_id}/${sample_id}_${suffix}"
    # create file as a flag to indicate successful completion

# fi
touch "ok.txt"



echo "\$(date)  COMPLETED!"
