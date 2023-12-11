#!/usr/bin/env bash
#SBATCH -J launch_apply_gene_expr_filter
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm_out/%x_%A.out
#SBATCH --error=slurm_out/%x_%A.out

job_min=1

# base_dir="${HOME}/Desktop/gaitigroup/Users"
base_dir="/cluster/projects/gaitigroup/Users"
work_dir=$base_dir/Joan/wes-mutation-calling


sample_dir="${work_dir}/output/test_set/mutations_prefiltered"
output_dir="${work_dir}/output/test_set/mutations_postprocessed"


# Determine job array limits
# A. Determine number of files
job_max=$(ls -d -R -- $sample_dir/*.txt | wc -l) 2>/dev/null
# B. Number of lines in a file
# job_max=$(wc -l < "${sample_ids}")
# job_max=1

echo $job_max

sbatch <<EOF
#!/usr/bin/env bash

#SBATCH -J postprocessing
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=joan.kant@uhn.ca
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=slurm_out/%x_%A_%a.out
#SBATCH --error=slurm_out/%x_%A_%a.out
#SBATCH --array=${job_min}-${job_max}

echo "Activating conda environment..."
source "\$HOME/miniforge3/bin/activate" "standard_env"

sample=\$(ls -d -- $sample_dir/*.txt | sed -n \${SLURM_ARRAY_TASK_ID}p)

Rscript "$work_dir/scripts/postprocessing.R" \
    --output_dir ${output_dir} \
    --input_file \${sample}
EOF
