project_dir="/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-mutation-calling"
wes_mutations="/Users/joankant/Library/CloudStorage/OneDrive-UHN/004_Projects/Lupus/NIHMS907788-supplement-10.xlsx"
rna_mutations_dir="${project_dir}/output/mutations"
output_dir="${project_dir}/output/figures"
sample_sheet="${project_dir}/misc/sample_sheet.csv"

source "/Users/joankant/miniforge-pypy3/bin/activate" "cci_env"


for rna_mutations in rna_mutations_dir/*;
do
echo $rna_mutations
Rscript "${project_dir}/scripts/comparison_wes_raw.R" \
    --wes_mutations ${wes_mutations} \
    --rna_mutations $rna_mutations \
    --output_dir ${output_dir} \
    --sample_sheet ${sample_sheet}
done
