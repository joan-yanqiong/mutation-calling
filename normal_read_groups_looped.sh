#!/bin/bash

#SBATCH -p all
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=4G 
#SBATCH --job-name read_groups_loop		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load picard

directory="/cluster/projects/gaitigroup/Users/Jahin/Results/bwa"
script_path="/cluster/projects/gaitigroup/Users/Jahin/Scripts/pipeline_scripts"

cd $directory
i=1
for run_id in SRR*
do
    i=`expr $i + 1`
    cd ${run_id}
    sbatch ${script_path}/normal_read_groups.sh $run_id $i
    cd ..
done