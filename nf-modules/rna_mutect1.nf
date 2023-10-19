process MUTECT_R1 {
    /*
    Run MuTect

    Input:
    ref_path: path to reference genome
    cosmic_vcf: path to COSMIC VCF
    dbSNP_vcf: path to dbSNP VCF
    run_id: Directory containing all current pipeline output files, incl. *_recal.bam

    */

    publishDir "${projectDir}/output/tumor", mode: "copy"

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

    publishDir "${projectDir}/output/tumor", mode: "copy"

    input:
    tuple val(ix),val(sample_id), path(mutect_vcf), path(mutect_vcf_index)

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_onco.maf"), path(mutect_vcf_index), emit: onco_maf
    path "${sample_id}/oncotator.log"

    script:
    template "oncotator.sh"

}