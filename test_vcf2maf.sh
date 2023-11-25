#!/usr/bin/env bash
# module load vep
# module load vcf2maf
module load samtools
vcf2maf_dir="/cluster/home/t119972uhn/downloads/mskcc-vcf2maf-754d68a"


# perl "/cluster/tools/software/centos7/vcf2maf/1.6.17/vcf2maf.pl"
input_vcf="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/test_folder/SRR5088818/SRR5088818_mutect_round1.hg19_multianno.vcf"
output_maf="/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/test_folder/SRR5088818/SRR5088818_mutect_round1.hg19_multianno.maf"

ref_fasta="/cluster/projects/gaitigroup/Users/Joan/001_data/Lupus/reference_genome/Homo_sapiens_assembly19.fasta"
# vep --offline --dir_cache $cache_dir
# ref_fasta="/cluster/tools/data/commondata/ensembl/vep/105/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
# vcf_filter="/cluster/tools/data/commondata/ensembl/vep/105/homo_sapiens/105_GRCh38/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz"

# cache_dir="/cluster/tools/data/commondata/ensembl/vep/105"
# vep_path="/cluster/tools/software/centos7/vep/105"

# .vep/homo_sapiens/95_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

# vep --offline --dir_cache $cache_dir --a GRCh38

# ref_fasta="/cluster/projects/gaitigroup/Users/Joan/001_data/Lupus/reference_genome/Homo_sapiens_assembly19.fasta"
perl ${vcf2maf_dir}/vcf2maf.pl \
    --input-vcf ${input_vcf} \
    --output-maf ${output_maf} \
    --inhibit-vep \
    --ref-fasta ${ref_fasta} \
    --tumor-id "SRR5088818"
    # --filter-vcf 0 \
    # --vep-data ${cache_dir} \
    # --vep-path ${vep_path} \
    # # --tumor-id "SRR5088818" \
    # # --inhibit-vef
