#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name mutect
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load mutect/1.1.5

interval_list=snp_mutations.intervals

mkdir -p ${tumor_id}
cd ${tumor_id}

java -jar \$mutect_dir/muTect.jar \
    -T MuTect \
    -R "../${ref_path}" \
    --input_file:tumor "../${tumor_hisat_bam}" \
    --input_file:normal "../${normal_tumor_hisat_bam}" \
    -cosmic "../${cosmic_vcf}" \
    -dbsnp "../${dbSNP_vcf}" \
    -o ${tumor_id}_second_call_stats.txt \
    -vcf ${tumor_id}_second_mutect.vcf \
    -U ALLOW_N_CIGAR_READS \
    --force_output \
    -L ${interval_list}
