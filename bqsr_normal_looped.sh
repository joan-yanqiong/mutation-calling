#!/bin/bash

#SBATCH -p all
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=4G
#SBATCH --job-name bqsr_single_loop		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"
directory=${prefix}/Results/bwa
ref_path=${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta
script_path="${prefix}/Scripts/pipeline_scripts"
dbSNP_vcf="${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.dbsnp138.vcf"

cd $directory

for run_id in SRR*
do
    cd $run_id
    sbatch ${script_path}/bqsr_normal.sh ${run_id} ${ref_path} ${dbSNP_vcf}
    cd ..
done