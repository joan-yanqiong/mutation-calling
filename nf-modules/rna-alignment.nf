process STAR_ALIGN {
    label "mem32"
    label "time_1h"
    /*
    Summary:

    Input:
    genome_dir - path to the directory (henceforth called ”genome directory”
    where the genome indices are stored
    fastq_dir - path to directory where the files containing the sequences to be
    mapped (e.g. RNA-seq FASTQ files) are stored

    Output:


    Ref:
    https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
    [p. 7-8]

    */
    publishDir "${projectDir}/${params.run_name}/output/tumor", mode: "symlink"

    input:
    path index_dir
    tuple val(ix), val(sample_id), path(fastq_dir)

    output:
    path "${sample_id}/mapped/${sample_id}_Log.final.out"
    path "${sample_id}/mapped/${sample_id}_Log.out"
    path "${sample_id}/mapped/${sample_id}__STARtmp"
    path "${sample_id}/mapped/${sample_id}_Log.progress.out"
    path "${sample_id}/mapped/${sample_id}_SJ.out.tab"
    tuple val(ix), val(sample_id), path("${sample_id}/mapped/${sample_id}_star_aligned.bam"), emit: output
    path "ok.txt"

    script:
    template "star_align.sh"
}
