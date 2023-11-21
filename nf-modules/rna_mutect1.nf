process MUTECT_R1 {
    label "time_8h"
    label "mem40"
    /*
    Run MuTect

    Input:
    ref_path: path to reference genome
    cosmic_vcf: path to COSMIC VCF
    dbSNP_vcf: path to dbSNP VCF
    run_id: Directory containing all current pipeline output files, incl. *_recal.bam

    */

    publishDir "${projectDir}/output/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(recal_bam), path(recal_bai)
    tuple path(ref_path), path(ref_path_dict), path(ref_path_fai)
    tuple path(dbSNP_vcf), path(dbSNP_vcf_idx)
    tuple path(cosmic_vcf), path(cosmic_vcf_idx)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_mutect.vcf"), path("${sample_id}/${sample_id}_mutect.vcf.idx"), emit: output
    path "${sample_id}/${sample_id}_call_stats.txt", emit: call_stats

    script:
    template "mutect_R1.sh"
}

process ONCOTATOR {
    label "time_30m"
    label "mem16"
    /*
    Summary: Annotate mutations obtained with Mutect

    Input:
    sample_id:
    mutect_vcf: {sample_id}_mutect.vcf
    mutect_vcf_idx: {sample_id}_mutect.vcf.idx

    Output:
    onco_maf: {sample_id}_onco.maf
    oncotator.log: oncotator.log

    Ref:
    */

    publishDir "${projectDir}/output/tumor", mode: "symlink"

    input:
    tuple val(ix),val(sample_id), path(mutect_vcf), path(mutect_vcf_index)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_onco.maf"), path(mutect_vcf_index), emit: onco_maf
    path "${sample_id}/oncotator.log"

    script:
    template "oncotator.sh"

}

// TODO: testing is necessary!
process CONVERT_TO_AVINPUT_R1 {
    label "time_30m"
    label "mem2"
    publishDir "${projectDir}/output/tumor", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(mutect_vcf), path(mutect_vcf)
    val suffix

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_${suffix}.avinput"), path(mutect_vcf_index)

    script:
    template "convert_to_avinput.sh"
}

process ANNOVAR_R1 {
    label "time_30m"
    label "mem2"
    publishDir "${projectDir}/output/tumor/", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(av_input), path(mutect_vcf_index)
    path human_db
    val suffix

    output:
    path "${sample_id}/${sample_id}_${suffix}.annovar.vcf"
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_${suffix}.annovar.hg19_multianno.txt"),
    path(mutect_vcf_index), emit: annovar

    script:
    template "annovar.sh"
}