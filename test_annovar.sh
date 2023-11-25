mutect_vcf="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_mutect.vcf"
human_db="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/data/humandb"
annovar="/cluster/tools/software/annovar/20180416"

sample_id="SRR5088818"
test_folder="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/test_folder"
suffix="mutect_round1"
module load annovar

echo "$(date)   Create output directory..."

echo "$(date)   Annotating variants..."

cd $test_folder

mkdir -p ${sample_id}

perl ${annovar}/table_annovar.pl \
                -vcfinput "${mutect_vcf}" \
                ${human_db} \
                -buildver hg19 \
                -remove \
                -protocol refGene \
                -operation g \
                -nastring . \
                -polish \
                --otherinfo \
                -out "${sample_id}/${sample_id}_${suffix}"
# for paralellization: --thread 4
echo "$(date)  COMPLETED!"


                # -out "${sample_id}/${sample_id}_${suffix}.annovar" \
