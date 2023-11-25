#!/usr/bin/env bash

#SBATCH -J test-funcotator
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
##SBATCH --partition=himem
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

echo "\$(date)   Loading modules..."
module load gatk

data_sources_path="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/001_data/funcotator_dataSources.v1.7.20200521s"
input_vcf="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_second_mutect.vcf"
sample_id="SRR5088818"
suffix="mutect_R2"
ref_path="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/001_data/reference_genome/Homo_sapiens_assembly19.fasta"
# echo "\$(date)  Creating output directory..."
# mkdir -p ${sample_id}
funcotator_path="/cluster/tools/software/centos7/gatk/4.4.0.0/gatk"

cd "/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/test_folder/SRR5088818"

echo "\$(date)  Use Funcotator for annotation..."
${funcotator_path} Funcotator \
    --variant ${input_vcf} \
    --reference ${ref_path} \
    --ref-version hg19 \
    --data-sources-path ${data_sources_path} \
    --output ${sample_id}_${suffix}.funcotated.maf \
    --output-file-format MAF

touch "ok.txt"

echo "\$(date)   COMPLETED!"