#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name mutect_loop
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load mutect/1.1.5

prefix="/cluster/projects/gaitigroup/Users/Jahin"
script_path="/cluster/projects/gaitigroup/Users/Jahin/Scripts/pipeline_scripts"

directory=${prefix}/Results/Practice/STAR
ref_path=${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta
# normal_path="${prefix}/Results/STAR/Normal/SRR5134767/reordered_normal.bam"
cosmic_vcf="${prefix}/Reference/b37_cosmic_v54_120711.vcf"
dbSNP_vcf="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.dbsnp138.vcf"

cd $directory

for run_id in SRR*
do
    cd ${run_id}
    sbatch ${script_path}/mutect.sh ${ref_path} ${cosmic_vcf} ${dbSNP_vcf} ${run_id}
    cd ..
done
