#!/bin/bash

#SBATCH -p himem
#SBATCH -t 0-10:00
#SBATCH -c 6
#SBATCH --mem=40G
#SBATCH --job-name apply_bqsr_loop		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"
directory=${prefix}/Results/bwa
ref_path=${prefix}/Reference/hg19-v0-Homo_sapiens_assembly19.fasta
script_path="${prefix}/Scripts/pipeline_scripts"

cd $directory

for run_id in SRR*
do
    cd $run_id
    sbatch ${script_path}/apply_bqsr_normal.sh ${run_id} ${ref_path}
    cd ..
done