process BWA_ALIGN {
    label "time_4h"
    label "mem12"
    /*
    Summary: Align normal sample to the reference genome using bwa

    Input:
    sample_id: Sample ID
    dir: Path to the directory containing the sample fastq files
    genome_dir: Path to the indexed reference genome

    Output:
    mapped_sam: Path to the mapped sam file

    Ref: https://bio-bwa.sourceforge.net/bwa.shtml
    */
    publishDir "${projectDir}/output/normal", mode: "symlink"

    input:
    tuple val(ix), val(sample_id), path(dir)
    path index_dir

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_mapped.sam"), emit: mapped_sam

    script:
    template "bwa_align.sh"
}
