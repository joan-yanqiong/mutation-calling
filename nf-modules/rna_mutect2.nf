process TUMOR_REALIGN_PREPROCESS {
    /*
    Summary:
    1. Generate intervals from maf or maflite
    2. Extracting paired reads from BAM
    3. Convert BAM to FASTQ

    Input:
    jar_files libdir
    sample_id: sample id, i.e. tumor_id
    maf: maf file that contains point mutation information, output from
    MUTECT(), e.g. {sample_id}_mutect.maf
    recal_bam: recalibrated bam file, output from APPLY_BQSR(), e.g. {sample_id}_recal.bam
    ref_path: path to reference genome (FASTA)

    Output:


    */
    label "time_8h"
    label "mem16"

    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"
    input:
    path jar_files
    tuple val(ix), val(sample_id), val(normal_id), val(pair_id), path(maf), path(maf_idx), path(recal_bam), path(recal_bai)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)
    val sample_type

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/snp_mutations.intervals"),
    path("${sample_id}/snp_mutations.intervals.bed"), emit: snp_mut_intervals
    path "${sample_id}/IDs_all.txt"
    path "${sample_id}/tmp_bam.bam"
    path "${sample_id}/tmp_header_T.sam"
    path "${sample_id}/tmp0_T.sam"
    path "${sample_id}/tmp_filteredbamT.sam"
    tuple val(ix), val(sample_id), val(normal_id), val(pair_id), path("${sample_id}/${sample_id}_tmp_sequence_1.fastq"), path("${sample_id}/${sample_id}_tmp_sequence_2.fastq"), val(sample_type), emit: output
    path "${sample_id}/ok.txt"

    shell:
    template "HiSat_realign_preprocess2.sh"
}


process NORMAL_REALIGN_PREPROCESS {
    label "time_8h"
    label "mem16"
    /*
    Input:
    jar_files libdir
    sample_id: sample id, should be of format {normal_id}_{tumor_id}
    maf: maf file that contains point mutation information, output from
    MUTECT(), e.g. {sample_id}_mutect.maf
    recal_bam: recalibrated bam file, output from APPLY_BQSR(), e.g. {sample_id}_recal.bam
    ref_path: path to reference genome (FASTA)

    Output:
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "symlink"

    input:
    path jar_files
    tuple val(ix), val(tumor_id), val(sample_id), val(pair_id), path(maf), path(maf_idx), path(recal_bam), path(recal_bai)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)
    val sample_type

    output:
    path "${sample_id}/snp_mutations.intervals"
    path "${sample_id}/snp_mutations.intervals.bed"
    path "${sample_id}/IDs_all.txt"
    path "${sample_id}/tmp_bam.bam"
    path "${sample_id}/tmp_header_T.sam"
    path "${sample_id}/tmp0_T.sam"
    path "${sample_id}/tmp_filteredbamT.sam"
    tuple val(ix), val(tumor_id), val(sample_id),  val(pair_id),
    path("${sample_id}/${pair_id}_tmp_sequence_1.fastq"),
    path("${sample_id}/${pair_id}_tmp_sequence_2.fastq"), val(sample_type), emit: output
    path "${sample_id}/ok.txt"


    shell:
    template "HiSat_realign_preprocess2.sh"
}

process TUMOR_HISAT_ALIGN {
    label "mem6"
    label "time_1h"
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    /*
    Summary: Align RNA-seq reads to reference genome using HISAT

    Input:
    hisat_index_dir: path to directory containing HISAT index files
    sample_id: sample id

    Output:
    hisat_sam: {sample_id}.aligned.sorted_by_coord.hisat2.sam

    Ref: http://daehwankimlab.github.io/hisat2/manual/
    */

    input:
    path index_dir
    tuple val(ix), val(tumor_id), val(normal_id), val(pair_id), path(fastq1), path(fastq2), val(sample_type)

    output:
    tuple val(ix), val(tumor_id), path("${tumor_id}/${tumor_id}_hisat2_aligned.sam"), emit: output
    path "ok.txt"

    script:
    template "hisat_align.sh"
}

process PAIR_HISAT_ALIGN {
    label "mem6"
    label "time_1h"
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "symlink"

    /*
    Summary: Align RNA-seq reads to reference genome using HISAT

    Input:
    hisat_index_dir: path to directory containing HISAT index files
    sample_id: {normal_id}_{tumor_id}

    Output:
    hisat_sam: {sample_id}.aligned.sorted_by_coord.hisat2.sam

    Ref: http://daehwankimlab.github.io/hisat2/manual/
    */

    input:
    path index_dir
    tuple val(ix), val(tumor_id), val(normal_id), val(pair_id), path(fastq1), path(fastq2), val(sample_type)

    output:
    tuple val(ix), val(normal_id), path("${normal_id}/${pair_id}_hisat2_aligned.sam"), emit: output
    path "ok.txt"

    script:
    template "hisat_align.sh"
}


process SORT_BAM_COORD_HIS_NORMAL{
    label "time_1h"
    label "mem6"
    /*
    Summary: Sorts the input SAM or BAM file by queryname

    Input:
    sample_id: Sample ID
    read_groups_bam: BAM file containing read groups

    Output:
    sample_id: Sample ID
    sorted_bam: Sorted BAM file

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/4418062801691-SortSam-Picard-
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(read_groups_sam)
    val sort_order

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_sortedBy_${sort_order}.bam"), emit: output
    path "ok.txt"

    script:
    template "sort_bam.sh"
}

process SORT_BAM_COORD_HIS_TUMOR{
    label "time_1h"
    label "mem6"

    /*
    Summary: Sorts the input SAM or BAM file by queryname

    Input:
    sample_id: Sample ID
    read_groups_bam: BAM file containing read groups

    Output:
    sample_id: Sample ID
    sorted_bam: Sorted BAM file

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/4418062801691-SortSam-Picard-
    */
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(read_groups_sam)
    val sort_order

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_sortedBy_${sort_order}.bam"), emit: output
    path "ok.txt"

    script:
    template "sort_bam.sh"
}

process INDEX_BAM_HIS_NORMAL {
    label "time_30m"
    label "mem2"
        /*
    Summary: Indexes a BAM file using samtools

    Input:
    sample_id: Sample ID
    sorted_bam: Sorted BAM file

    Output:
    sample_id: Sample ID
    sorted_bam: Sorted BAM file

    Ref: http://www.htslib.org/doc/samtools-index.html
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(bam_file)

    output:
    tuple val(ix), val(sample_id), path(bam_file), path("${sample_id}/${bam_file}.bai"), emit: output
    path "ok.txt"

    script:
    template "index_bam.sh"
}


process INDEX_BAM_HIS_TUMOR {
    label "time_30m"
    label "mem2"
    /*
    Summary: Indexes a BAM file using samtools

    Input:
    sample_id: Sample ID
    sorted_bam: Sorted BAM file

    Output:
    sample_id: Sample ID
    sorted_bam: Sorted BAM file

    Ref: http://www.htslib.org/doc/samtools-index.html
    */
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(bam_file)

    output:
    tuple val(ix), val(sample_id), path(bam_file), path("${sample_id}/${bam_file}.bai"), emit: output
    path "ok.txt"

    script:
    template "index_bam.sh"
}
process MUTECT_R2 {
    label "time_8h"
    label "mem40"
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    /*
    Summary: Run MuTect2 on tumor-normal pair

    Input:
    ref_path: path to reference genome (FASTA), e.g. hg19-v0-Homo_sapiens_assembly19.fasta
    ref_gtf: path to reference GTF file, e.g. Homo_sapiens_assembly19.gtf
    cosmic_vcf: path to COSMIC VCF file, e.g. b37_cosmic_v54_120711.vcf
    dbsnp_vcf: path to dbSNP VCF file, e.g.
    hg19-v0-Homo_sapiens_assembly19.dbsnp138.vcf

    tumor_id
    normal_id

    tumor_hisat_bam: ${tumor_id}.aligned.sorted_by_coord.hisat2.bam
    normal_tumor_hisat_bam:
    ${normal_id}/${normal_id}_${tumor_id}.aligned.sorted_by_coord.hisat2.bam

    Output:

    */
    input:
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)
    tuple path(dbSNP_vcf), path(dbSNP_vcf_idx)
    tuple path(cosmic_vcf), path(cosmic_vcf_idx)
    tuple val(ix), val(sample_id), val(_normal_id), val(_pair_id), path(tumor_hisat_bam), path(tumor_hisat_bam_bai), path(normal_tumor_hisat_bam), path(normal_tumor_hisat_bam_bai), path(snp_mut_intervals), path(snp_mut_intervals_bed)
    val suffix

    output:
    path "${sample_id}/${sample_id}_${suffix}_call_stats.txt"
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_${suffix}.vcf"), path("${sample_id}/${sample_id}_${suffix}.vcf.idx"), emit: output
    path "ok.txt"

    script:
    template "mutect_R2.sh"
}

process FILTERING {
    label "time_30m"
    label "mem2"
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"
    /*
    SUMMARY:
    Filtering of results from variant calling (Mutect Round 2)
    1. Keep only sites that are PASS in Mutect
    2. Only keep sites that have a min. number of reads supporting these
       alternate allele
    3. Remove sites that are in ExAC, DARNED, RADAR, pseudo-genes

    Input:
    min_alt_count: the minimal number of reads required supporting the alternate allele
    exac_bed: originally downloaded on June 2015 https://doi.org/10.1038/nature19057
    darned_bed: originally downloaded on September 2015 https://doi.org/10.1093/nar/gks961
    radar_bed: originally downloaded on September 2015 https://doi.org/10.1093/nar/gkt996

    Additional information:
    exac_bed, darned_bed and radar_bed downloaded from  https://github.com/kyizhak/RNA_MUTECT_1.0
    and converted to BED format using Python/convert_to_bed.ipynb

    pseudo_genes_bed downloaded from
    http://pseudofam.pseudogene.org/pages/info/download.jsf

    pseudo_genes_bed2 constructed using bioMart, see scripts/get_annotations.R
    */
    input:
    tuple val(ix), val(sample_id), path(mutect_vcf), path(mutect_vcf_idx)
    path exac_bed
    path darned_bed
    path radar_bed
    path pseudo_genes_bed
    path pseudo_genes_bed2
    val min_alt_counts

    output:
    path("${sample_id}/mutations_filtered/${sample_id}.pass.recode.vcf")
    path("${sample_id}/mutations_filtered/${sample_id}.mac.recode.vcf")
    path("${sample_id}/mutations_filtered/${sample_id}.exac.recode.vcf")
    path("${sample_id}/mutations_filtered/${sample_id}.darned.recode.vcf")
    path("${sample_id}/mutations_filtered/${sample_id}.radar.recode.vcf")
    path("${sample_id}/mutations_filtered/${sample_id}.a_pseudo_genes.recode.vcf")
    tuple val(ix), val(sample_id), path("${sample_id}/mutations_filtered/${sample_id}.b_pseudo_genes.recode.vcf"), emit: output
    path "${sample_id}/mutations_filtered/ok.txt"

    script:
    template "filtering.sh"
}

process ANNOVAR_R2 {
    label "time_30m"
    label "mem2"
    publishDir "${projectDir}/output/${params.run_name}/tumor/", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(mutct_vcf), path(mutect_vcf_index)
    path human_db
    val suffix

    output:
    path "${sample_id}/${sample_id}_${suffix}.avinput"
    path "${sample_id}/${sample_id}_${suffix}.hg19_multianno.txt"
    path "${sample_id}/${sample_id}_${suffix}.hg19_multianno.vcf"
    path "ok.txt"

    script:
    template "annovar.sh"
}

process FUNCOTATOR {
    label "mem16"
    label "time_30m"
    publishDir "${projectDir}/output/${params.run_name}/mutations_prefiltered", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(input_vcf)
    path data_source_path
    val suffix
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    path "${sample_id}/${sample_id}_${suffix}.funcotated.maf"
    path "ok.txt"

    script:
    template "funcotator.sh"
}