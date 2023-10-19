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
    label "realign_preprocess"

    publishDir "${projectDir}/output/tumor", mode: "copy"

    input:
    path jar_files
    tuple val(ix), val(sample_id), val(normal_id), val(pair_id), path(maf), path(maf_idx), path(recal_bam), path(recal_bai)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/snp_mutations.intervals"),
    path("${sample_id}/snp_mutations.intervals.bed"), emit: snp_mut_intervals
    path "${sample_id}/IDs_all.txt"
    path "${sample_id}/tmp_bam.bam"
    path "${sample_id}/tmp_header_T.sam"
    path "${sample_id}/tmp0_T.sam"
    path "${sample_id}/tmp_filteredbamT.sam"
    tuple val(sample_id), val(normal_id), val(pair_id), path("${sample_id}/${sample_id}_tmp_sequence_1.fastq"), path("${sample_id}/${sample_id}_tmp_sequence_2.fastq"), val("tumor"), emit: output

    shell:
    template "HiSat_realign_preprocess2.sh"
}


process NORMAL_REALIGN_PREPROCESS {
    label "realign_preprocess"

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
    publishDir "${projectDir}/output/normal", mode: "copy"

    input:
    path jar_files
    tuple val(ix), val(tumor_id), val(normal_id), val(pair_id), path(maf), path(maf_idx), path(recal_bam), path(recal_bai)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    path "${normal_id}/snp_mutations.intervals"
    path "${normal_id}/snp_mutations.intervals.bed"
    path "${normal_id}/IDs_all.txt"
    path "${normal_id}/tmp_bam.bam"
    path "${normal_id}/tmp_header_T.sam"
    path "${normal_id}/tmp0_T.sam"
    path "${normal_id}/tmp_filteredbamT.sam"
    tuple val(ix), val(tumor_id), val(normal_id),  val(pair_id),
    path("${normal_id}/${pair_id}_tmp_sequence_1.fastq"),
    path("${normal_id}/${pair_id}_tmp_sequence_2.fastq"), val("normal"), emit: output

    shell:
    template "HiSat_realign_preprocess2.sh"
}

process TUMOR_REALIGN_ANNOTATE {
    publishDir "${projectDir}/output/tumor", mode: "copy"

    label "realign_annotate"
    /*
    Input:
    sample_id: tumor_id
    */
    input:
    tuple val(sample_id), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path(fastq1), path(fastq2), emit:output
    path "${sample_id}.rna_reads_fastq_list.list"

    shell:
    """
    #!/bin/bash

    mkdir -p !{sample_id}
    cd !{sample_id}

    echo "Writing out annotation file for upload"
    echo -e "!{fastq1}\n!{fastq2}"   >> !{sample_id}.rna_reads_fastq_list.list
    echo "COMPLETED!"
    """
}

process NORMAL_REALIGN_ANNOTATE {
    label "realign_annotate"

    publishDir "${projectDir}/output/normal", mode: "copy"

    input:
    tuple val(normal_id), val(tumor_id), val(pair_id), path(fastq1), path(fastq2)

    output:
    tuple val(normal_id), val(tumor_id), val(pair_id), path(fastq1), path(fastq2), emit: output
    path "${pair_id}.rna_reads_fastq_list.list"

    shell:
    """
    #!/bin/bash

    mkdir -p !{normal_id}
    cd !{normal_id}

    echo "writing out annotation file for upload"
    echo -e "!{fastq1}\n!{fastq2}"   >> !{pair_id}.rna_reads_fastq_list.list
    echo "COMPLETED!"

    """

}

process TUMOR_HISAT_ALIGN {
    label "hisat_align"
    publishDir "${projectDir}/output/tumor", mode: "copy"

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
    tuple val(ix), val(tumor_id), path("${tumor_id}/${tumor_id}.aligned.sorted_by_coord.hisat2.sam")

    script:
    template "hisat_align.sh"
}

process PAIR_HISAT_ALIGN {
    label "hisat_align"
    publishDir "${projectDir}/output/normal", mode: "copy"

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
    tuple val(ix), path(tumor_id), path("${normal_id}/${pair_id}.aligned.sorted_by_coord.hisat2.sam")

    script:
    template "hisat_align.sh"
}

process MUTECT_R2 {
    publishDir "${projectDir}/output/tumor", mode: "copy"

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
    tuple val(ix), val(sample_id), path(tumor_hisat_bam), path(normal_tumor_hisat_bam), path(snp_mut_intervals), path(snp_mut_intervals_bed)

    output:
    path "${sample_id}/${sample_id}_second_call_stats.*"
    path "${sample_id}/${sample_id}_second_mutect.vcf.*"

    script:
    template "mutect_R2.sh"
}