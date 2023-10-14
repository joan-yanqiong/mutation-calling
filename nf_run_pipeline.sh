#!/usr/bin/env bash
#SBATCH -J run_nf_mutation_calling
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=long
##SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=21-00:00:00
##SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

module load java
# module load apptainer

# Local
# base_dir="/Users/joankant/Desktop/gaitigroup/Users/Joan"
# nf_exec="/Users/joankant/Library/CloudStorage/OneDrive-UHN/nextflow"
# work_dir="${base_dir}/nf_work"
# nf_profile="standard"

# H4H
base_dir="/cluster/projects/gaitigroup/Users/Joan/"
nf_exec="/cluster/home/t119972uhn/nextflow"
work_dir="${base_dir}/nf_work"
nf_profile="slurm"

# Create work directory if not existing
mkdir -p $work_dir

project_dir="${base_dir}/h4h-mutation-calling"

# USER PARAMETERS
sample_sheet="${project_dir}/misc/sample_sheet_test.csv"

# HELPER FILES (REFERENCES)
ref_genome="${project_dir}/data/reference_genome/Homo_sapiens_assembly19.fasta"
gtf_path="${project_dir}/data/References/Homo_sapiens_assembly19.gtf"

dbSNP="${project_dir}/data/References/hg19-v0-Homo_sapiens_assembly19.dbsnp138.vcf"
dbcosmic="${project_dir}/data/References/b37_cosmic_v54_120711.vcf"

indel_db1="${project_dir}/data/References/1000G_phase1.indels.b37.vcf"
indel_db2="${project_dir}/data/References/Mills_and_1000G_gold_standard.indels.b37.vcf"

jar_files="${project_dir}/data/jar_files"

# Create main output directories
mkdir -p "${project_dir}/data/indices"
mkdir -p "${project_dir}/logs"

# Start the pipeline
${nf_exec} run ${project_dir} \
    -resume "c19daa3f-d1be-4047-ab45-f90f595f53f5" -profile ${nf_profile} \
    -w $work_dir \
    --sample_sheet ${sample_sheet} \
    --ref_genome ${ref_genome} \
    --gtf_path ${gtf_path} \
    --dbSNP ${dbSNP} \
    --dbcosmic ${dbcosmic} \
    --indel_db1 ${indel_db1} \
    --indel_db2 ${indel_db2} \
    --jar_files ${jar_files} \
    --outdir ${project_dir}/logs



# /cluster/projects/gaitigroup/Users/Jahin/Results/Practice/STAR/SRR5088818/snp_mutations.intervals /cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/tmp/snp_mutations.intervals

# /cluster/projects/gaitigroup/Users/Jahin/Results/Practice/STAR/SRR5088818/snp_mutations.intervals.bed /cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/tmp/snp_mutations.intervals.bed