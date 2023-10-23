echo "\$(date)  Loading mutect module..."
module load mutect/1.1.5

echo "\$(date)  Creating directory for output..."
mkdir -p ${sample_id}

echo "\$(date) Enter created folder..."
cd ${sample_id}

echo "\$(date)  Run mutect..."
java -jar \$mutect_dir/muTect.jar \
    -T MuTect \
    -R "../${ref_path}" \
    --input_file:tumor "../${tumor_hisat_bam}" \
    --input_file:normal "../${normal_tumor_hisat_bam}" \
    -cosmic "../${cosmic_vcf}" \
    -dbsnp "../${dbSNP_vcf}" \
    -o ${sample_id}_second_call_stats.txt \
    -vcf ${sample_id}_second_mutect.vcf \
    -U ALLOW_N_CIGAR_READS \
    --force_output \
    -L "../${snp_mut_intervals}"

echo "\$(date)  COMPLETED!"
