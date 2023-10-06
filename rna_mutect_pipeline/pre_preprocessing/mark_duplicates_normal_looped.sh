#!/bin/bash

#SBATCH -p all
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=2G 
#SBATCH --job-name mark_duplicates_and_split_cigars_loop		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

prefix="/cluster/projects/gaitigroup/Users/Jahin"
directory=${prefix}/Results/bwa
script_path="/cluster/projects/gaitigroup/Users/Jahin/Scripts/pipeline_scripts"

cd $directory

for run_id in SRR*
do
    cd ${run_id}
    sbatch ${script_path}/mark_duplicates_normal.sh $run_id
    cd ..
done