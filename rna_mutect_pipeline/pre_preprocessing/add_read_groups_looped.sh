#!/bin/bash

#SBATCH -p all
#SBATCH -t 0-6:00
#SBATCH -c 6
#SBATCH --mem=4G 
#SBATCH --job-name read_groups_loop		
#SBATCH -o /cluster/projects/gaitigroup/Users/Jahin/Logs/%x-%j.out

module load picard/2.10.9

directory=$1 # output directory
# csv_path="/cluster/projects/gaitigroup/Users/Jahin/SRR_to_GSM_mappings.csv" 
script_path="/cluster/projects/gaitigroup/Users/Jahin/Scripts/pipeline_scripts"

cd $directory
i=1
for run_id in SRR*
do
    cd ${run_id}
    sbatch ${script_path}/add_read_groups.sh $i $run_id
    i=`expr $i + 1`
    cd ..
done
