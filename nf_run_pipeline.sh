#!/usr/bin/env bash
#SBATCH -J run_nf_mutation_calling
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=21-00:00:00
#SBATCH --output=slurm_out/%x_%j.out
#SBATCH --error=slurm_out/%x_%j.out

module load java

# H4H
base_dir="/cluster/projects/gaitigroup/Users/Joan/"
nf_exec="/cluster/home/t119972uhn/nextflow"
work_dir="${base_dir}/nf_work"
nf_profile="slurm"

echo "$(date)   Create work directory for nextflow if not existing..."
mkdir -p "${work_dir}"

echo "$(date)   Setup paths..."
project_dir="${base_dir}/h4h-mutation-calling"
sample_sheet="${project_dir}/misc/sample_sheet_subset.csv"

# HELPER FILES (REFERENCES)
ref_genome="${project_dir}/data/reference_genome/Homo_sapiens_assembly19.fasta"
gtf_path="${project_dir}/data/References/Homo_sapiens_assembly19.gtf"

# Databases
dbcosmic="${project_dir}/data/References/b37_cosmic_v54_120711.vcf"
dbSNP="${project_dir}/data/References_v2/Homo_sapiens_assembly19.dbsnp138.vcf"
indel_db1="${project_dir}/data/References_v2/1000G_omni2.5.b37.vcf.gz"
indel_db2="${project_dir}/data/References_v2/1000G_phase1.snps.high_confidence.b37.vcf.gz"

human_db="${project_dir}/data/humandb/"

# Databases for removing gene sites
exac_bed="${project_dir}/data/sites_to_remove_dbs/exac_mat.bed"
darned_bed="${project_dir}/data/sites_to_remove_dbs/Darned_mat.bed"
radar_bed="${project_dir}/data/sites_to_remove_dbs/Radar_mat.bed"
pseudo_genes_bed="${project_dir}/data/sites_to_remove_dbs/pseudo_genes.bed"
pseudo_genes_bed2="${project_dir}/data/sites_to_remove_dbs/pseudo_genes_2.bed"

# Required tools
jar_files="${project_dir}/data/jar_files"

# Indexing with star, hisat and bwa
star_index_dir="${project_dir}/data/indices/STAR_Homo_sapiens_assembly19"
hisat_index_dir="/cluster/projects/gaitigroup/Users/Joan/utilities/skeleton/nf_pipeline/data/indices/HISAT_Homo_sapiens_assembly19"
bwa_index_dir="${project_dir}/data/indices/BWA_Homo_sapiens_assembly19"

# Create main output directories
echo "$(date)   Create output directories..."
mkdir -p "${project_dir}/data/indices"
mkdir -p "${project_dir}/logs"
mkdir -p "${project_dir}/output/normal"
mkdir -p "${project_dir}/output/tumor"
mkdir -p "${project_dir}/output/mutations_prefiltered"

# Start the pipeline
echo "$(date)   Start the pipeline..."
${nf_exec} run ${project_dir} -with-report -with-trace \
    -profile ${nf_profile} \
    -resume "665f4183-3717-4b49-9089-a99ba7b10b0f" \
    -w $work_dir \
    --sample_sheet ${sample_sheet} \
    --ref_genome ${ref_genome} \
    --gtf_path ${gtf_path} \
    --dbSNP ${dbSNP} \
    --dbcosmic ${dbcosmic} \
    --indel_db1 ${indel_db1} \
    --indel_db2 ${indel_db2} \
    --jar_files ${jar_files} \
    --outdir ${project_dir}/logs \
    --star_index_dir ${star_index_dir} \
    --hisat_index_dir ${hisat_index_dir} \
    --bwa_index_dir ${bwa_index_dir} \
    --human_db ${human_db} \
    --exac_bed ${exac_bed} \
    --darned_bed ${darned_bed} \
    --radar_bed ${radar_bed} \
    --pseudo_genes_bed ${pseudo_genes_bed} \
    --pseudo_genes_bed2 ${pseudo_genes_bed2} \
    --min_alt_counts ${min_alt_counts}

echo "$(date)   COMPLETED!"