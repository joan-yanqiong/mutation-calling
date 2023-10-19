#!/usr/bin/env bash
#SBATCH -J run_nf_mutation_calling
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=long
##SBATCH --partition=all
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=21-00:00:00
##SBATCH --time=00:10:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

base_dir="/cluster/projects/gaitigroup/Users/Joan/"
nf_exec="/cluster/home/t119972uhn/nextflow"
work_dir="${base_dir}/nf_work"
nf_profile="slurm"

# HELPER FILES (REFERENCES)
ref_genome="${project_dir}/data/reference_genome/Homo_sapiens_assembly19.fasta"
gtf_path="${project_dir}/data/References/Homo_sapiens_assembly19.gtf"

dbSNP="${project_dir}/data/References/hg19-v0-Homo_sapiens_assembly19.dbsnp138.vcf"
dbcosmic="${project_dir}/data/References/b37_cosmic_v54_120711.vcf"

indel_db1="${project_dir}/data/References/1000G_phase1.indels.b37.vcf"
indel_db2="${project_dir}/data/References/Mills_and_1000G_gold_standard.indels.b37.vcf"

jar_files="${project_dir}/data/jar_files"


bam_file="/cluster/projects/gaitigroup/Users/Jahin/Results/bwa/SRR5134767/SRR5134767_read_groups.bam"

bash "templates/sort_bam.sh"