process ADD_READ_GROUPS {
    label 'mem2'
    label 'time_1h'
    /* ADD READ GROUPS
    */

    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(mapped_bam)
    val(sample_type)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_RG.bam"), emit: output
    path "ok.txt"

    script:
    template "add_read_groups.sh"
}

process MARK_DUPLICATES {
    label 'mem4'
    label 'time_2h'
    /*
    Summary: A better duplication marking algorithm that handles all cases including clipped and gapped alignments.

    Input:
    ix, sample_id, dir: sample information from sample sheet
    input_bam: Directory containing a BAM/SAM/CRAM file containing reads, here _read_groups

    Output:
    - marked_dup_bam: The output file to write marked records to
    - metrics_file: File to write duplication metrics to

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
    */
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(read_groups_file)

    output:
    path "${sample_id}/${sample_id}_marked_dup_metrics.txt", emit: metrics_file
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_marked_dup.bam"), emit: output
    path "ok.txt"

    script:
    template "mark_duplicates.sh"
}

process SPLIT_CIGARS {
    label 'mem8'
    label 'time_4h'
    /*
    Input:
    ix, sample_id, dir: information from sample sheet
    marked_dup_bam: Directory containing a BAM/SAM/CRAM file containing reads.
    ref_path: Reference sequence file

    Output:
    split_bam: BAM file with reads split at N CIGAR elements and
    CIGAR strings updated.
    split_bai

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads

    */
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(marked_dup_bam)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_split.bam"),
    path("${sample_id}/${sample_id}_split.bai"), emit: output
    path "${sample_id}/ok.txt"

    script:
    template "split_cigars.sh"
}


process BQSR_TABLE {
    label "time_8h"
    label "mem8"
    /*
    Summary: Generates recalibration table for Base Quality Score Recalibration (BQSR)

    Input:
    ix, sample_id, dir: information from sample sheet
    bam_file: {sample_id}_split.bam
    ref_path: Reference sequence file (FASTA), incl. fai and dict
    dbSNP_vcf: dbSNP VCF file, incl index

    Output:
    bqsr_table: The output recalibration table file to create

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
    */
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(bam_file), path(bai_file)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)
    tuple path(dbSNP_vcf), path(dbSNP_vcf_idx)

    output:
    tuple val(ix), val(sample_id), path(bam_file), path("${sample_id}/${sample_id}_recal_data.table"), emit: output
    path "${sample_id}/ok.txt"

    script:
    template "bqsr_table.sh"
}

process APPLY_BQSR {
    label "time_8h"
    label "mem8"
    /*
    Summary: Apply base quality score recalibration

    Input:
    bam_file: {sample_id}_split.bam
    bqsr_table: {sample_id}_recal_data.table
    ref_path: Reference sequence file (FASTA)

    Output:
    recal_data_table: A BAM or CRAM file containing the recalibrated read data

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR

    */
    publishDir "${projectDir}/output/${params.run_name}/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(bam_file), path(recal_data_table)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_recal.bam"), path("${sample_id}/${sample_id}_recal.bai"), emit: output
    path "${sample_id}/ok.txt"


    script:
    template "apply_bqsr.sh"
}
