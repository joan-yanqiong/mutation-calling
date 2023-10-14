process HISAT_ALIGN_TUMOR_BAM {
    label "sam_to_bam"

    /*
    Summary: Convert SAM to BAM file

    Input:
    sample_id: Sample ID
    sam_file: SAM file containing reads

    Output:
    sample_id: Sample ID
    bam_file: BAM file containing reads

    */
    publishDir "${projectDir}/output/tumor", mode: "copy"

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}/${sam_file.simpleName}.bam")

    script:
    template "sam_to_bam.sh"
}

process HISAT_ALIGN_PAIR_BAM {
    label "sam_to_bam"
    /*
    Summary: Convert SAM to BAM file

    Input:
    sample_id: Sample ID
    sam_file: SAM file containing reads

    Output:
    sample_id: Sample ID
    bam_file: BAM file containing reads

    */
    publishDir "${projectDir}/output/normal", mode: "copy"

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}/${sam_file.simpleName}.bam")

    script:
    template "sam_to_bam.sh"
}