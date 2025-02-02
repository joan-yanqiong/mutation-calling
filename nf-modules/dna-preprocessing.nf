process ADD_READ_GROUPS_NORMAL {
    /* ADD READ GROUPS
    */
    label 'mem2'
    label 'time_30m'

    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(mapped_bam)
    val sample_type


    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_RG.bam"), emit: output
    path "ok.txt"

    script:
    template "add_read_groups.sh"
}

process SORT_BAM {
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
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(read_groups_sam)
    val sort_order

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_sortedBy_${sort_order}.bam"), emit: output
    path "ok.txt"

    script:
    template "sort_bam.sh"
}

process MARK_DUPLICATES_NORMAL {
    label 'mem4'
    label 'time_2h'

    /*
    Summary: A better duplication marking algorithm that handles all cases including clipped and gapped alignments.

    Input:
    ix, sample_id, dir: sample information from sample sheet
    input_bam: Directory containing a BAM/SAM/CRAM file containing reads, here
    sorted bam file

    Output:
    - marked_dup_bam: The output file to write marked records to
    - metrics_file: File to write duplication metrics to

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(read_groups_file)

    output:
    path "${sample_id}/${sample_id}_marked_dup_metrics.txt", emit: metrics_file
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_marked_dup.bam"), emit: output
    path "ok.txt"

    script:
    template "mark_duplicates.sh"
}


process SORT_BAM_COORD {
    label "time_1h"
    label "mem4"
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
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(read_groups_sam)
    val sort_order

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_sortedBy_${sort_order}.bam"), emit: output
    path "ok.txt"

    script:
    template "sort_bam.sh"
}

process INDEX_BAM {
    label "time_10m"
    label "mem1"
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
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(bam_file)

    output:
    tuple val(ix), val(sample_id), path(bam_file), path("${sample_id}/${bam_file.name}.bai"), emit: output
    path "ok.txt"

    script:
    template "index_bam.sh"
}

process INDEL_REALIGN_TARGET {
    label "time_30m"
    label "mem6"
    /*
    Summary:
    1. Sort by coord,
    2. index bam (creates a BAI index)
    3. create targets for realigning around indels

    Input:
    sample_id: Sample ID
    marked_dup_bam: BAM file containing marked duplicates
    ref_path: Path to reference genome
    indel_db1: Path to first indel database
    indel_db2: Path to second indel database

    Output

    Ref:
    https://gatk.broadinstitute.org/hc/en-us/articles/360036510732-SortSam-Picard-
    http://www.htslib.org/doc/samtools-index.html
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(marked_dup_bam), path(marked_dup_bam_bai)
    tuple path(indel_db1), path(indel_db1_idx)
    tuple path(indel_db2), path(indel_db2_idx)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    tuple val(ix), val(sample_id), path(marked_dup_bam), path(marked_dup_bam_bai), path("${sample_id}/${sample_id}_realigner.intervals"), emit: output
    path "ok.txt"

    script:
    template "indel_realign_target.sh"

}

process INDEL_REALIGNER {
    label "time_2h"
    label "mem4"
    /*
    Summary: Realign indels using targets from previous step

    Input:
    sample_id: Sample ID
    marked_dup_bam: BAM file containing marked duplicates sorted by coordinates
    realigner_intervals: File containing intervals to realign around
    ref_path: Path to reference genome
    indel_db1: Path to first indel database
    indel_db2: Path to second indel database

    Ref: https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Perform_local_realignment_around_indels.md#section2
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(marked_dup_bam), path(marked_dup_bai), path(realigner_intervals)
    tuple path(indel_db1), path(indel_db1_idx)
    tuple path(indel_db2), path(indel_db2_idx)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    tuple val(ix), val(sample_id),
    path("${sample_id}/${sample_id}_realigned.bam"), path("${sample_id}/${sample_id}_realigned.bai"), emit: output
    path "ok.txt"

    script:
    template "indel_realigner.sh"

}

process NORMAL_BQSR_TABLE {
    label "time_1h"
    label "mem2"
    /*
    Summary: Generates recalibration table for Base Quality Score Recalibration (BQSR)

    Input:
    ix, sample_id, dir: information from sample sheet
    bam_file: {sample_id}_realigned.bam
    ref_path: Reference sequence file (FASTA)
    dbSNP_vcf: dbSNP VCF file

    Output:
    bqsr_table: The output recalibration table file to create

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(bam_file), path(bai_file)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)
    tuple path(dbSNP_vcf), path(dbSNP_vcf_idx)

    output:
    tuple val(ix), val(sample_id), path(bam_file), path("${sample_id}/${sample_id}_recal_data.table"), emit: output
    path "ok.txt"

    script:
    template "bqsr_table.sh"
}

process NORMAL_APPLY_BQSR {
    label "time_1h"
    label "mem1"
    /*
    Summary: Apply base quality score recalibration

    Input:
    bam_file: {sample_id}_realigned.bam
    bqsr_table: {sample_id}_recal_data.table
    ref_path: Reference sequence file (FASTA)

    Output:
    recal_data_table: A BAM or CRAM file containing the recalibrated read data

    Ref: https://gatk.broadinstitute.org/hc/en-us/articles/360037055712-ApplyBQSR
    */
    publishDir "${projectDir}/output/${params.run_name}/normal", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(bam_file), path(recal_data_table)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_recal.bam"), path("${sample_id}/${sample_id}_recal.bai"), emit: output
    path "${sample_id}/ok.txt"

    script:
    template "apply_bqsr.sh"
}