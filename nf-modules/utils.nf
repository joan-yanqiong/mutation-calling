process HISAT_ALIGN_TUMOR_BAM {
    label "time_1h"
    label "mem8"
    /*
    Summary: Convert SAM to BAM file

    Input:
    sample_id: Sample ID
    sam_file: SAM file containing reads

    Output:
    sample_id: Sample ID
    bam_file: BAM file containing reads

    */
    publishDir "${projectDir}/${params.run_name}/output/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(sam_file)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sam_file.simpleName}.bam")

    script:
    template "sam_to_bam.sh"
}

process HISAT_ALIGN_PAIR_BAM {
    label "time_1h"
    label "mem8"
    /*
    Summary: Convert SAM to BAM file

    Input:
    sample_id: Sample ID
    sam_file: SAM file containing reads

    Output:
    sample_id: Sample ID
    bam_file: BAM file containing reads

    */
    publishDir "${projectDir}/${params.run_name}/output/normal", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(sam_file)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sam_file.simpleName}.bam")

    script:
    template "sam_to_bam.sh"
}

process SAM_TO_BAM {
    label "time_1h"
    label "mem8"
    /*
    Summary: Convert SAM to BAM file

    Input:
    sample_id: Sample ID
    sam_file: SAM file containing reads

    Output:
    sample_id: Sample ID
    bam_file: BAM file containing reads

    */
    publishDir "${projectDir}/${params.run_name}/output/normal", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(sam_file)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_mapped.bam")

    script:
    template "sam_to_bam_sq.sh"
}

// process CONVERT_TO_MAF {
//     label "very_short_process"
//     publishDir "${projectDir}/${params.run_name}/output/tumor", mode: "symlink"

//     input:
//     tuple val(ix), val(sample_id), path(input_file)
//     val suffix

//     output:
//     tuple val(ix), val(sample_id),  path("${sample_id}/${sample_id}_${suffix}.maf")

//     script:
//     """
//     Rscript "${projectDir}/scripts/annovar_to_maf.R" \
//         --input_file "\$PWD/${input_file}" \
//         --output_dir "\$PWD/${sample_id}" \
//         --output_file "\$PWD/${sample_id}/${sample_id}_${suffix}"
//     """

// }