#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name bqsr_table_loop
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"
directory=${prefix}/Results/Practice/STAR
ref_path=${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta
script_path="${prefix}/Scripts/pipeline_scripts"
dbSNP_vcf="${prefix}/Reference/common_all_20180423.vcf" # this needs to change

cd $directory

for run_id in SRR*
do
    cd $run_id
    sbatch ${script_path}/bqsr_table.sh ${ref_path} ${dbSNP_vcf} ${run_id}
    cd ..
done
