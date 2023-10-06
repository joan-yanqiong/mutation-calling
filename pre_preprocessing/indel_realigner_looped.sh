#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name indel_realigner_loop	
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"
directory=${prefix}/Results/bwa
ref_path=${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta
indel_db1="${prefix}/Reference/1000G_phase1.indels.b37.vcf"
indel_db2="${prefix}/Reference/Mills_and_1000G_gold_standard.indels.b37.vcf"
script_path="${prefix}/Scripts/pipeline_scripts"

cd $directory

for run_id in SRR*
do
    cd $run_id
    sbatch ${script_path}/indel_realigner.sh ${run_id} ${ref_path} ${indel_db1} ${indel_db2}
    cd ..
done