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

    shell:
    template "HiSat_realign_preprocess2.sh"
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
    tuple val(ix), val(tumor_id), path("${tumor_id}/${tumor_id}_aligned_hisat2.sam")

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
    tuple val(ix), val(normal_id), path("${normal_id}/${pair_id}_aligned_hisat2.sam")

    script:
    template "hisat_align.sh"
}


process SORT_BAM_COORD_HIS_NORMAL{
    label "very_short_process"

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
    publishDir "${projectDir}/output/tumor", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(read_groups_sam)
    val sort_order
    val suffix

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_${suffix}.bam"), emit: output

    script:
    template "sort_bam.sh"
}

process SORT_BAM_COORD_HIS_TUMOR{
    label "very_short_process"

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
    publishDir "${projectDir}/output/tumor", mode: "copy"

    input:
    tuple val(ix), val(sample_id), path(read_groups_sam)
    val sort_order
    val suffix

    output:
    tuple val(ix), val(sample_id), path("${sample_id}/${sample_id}_${suffix}.bam"), emit: output

    script:
    template "sort_bam.sh"
}

process INDEX_BAM_HIS_NORMAL {
    label "very_short_process"
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
    publishDir "${projectDir}/output/normal", mode: "copy"

    input:
   tuple val(ix), val(sample_id), path(bam_file)

    output:
    tuple val(ix), val(sample_id), path(bam_file), path("${sample_id}/${bam_file}.bai"), emit: output

    script:
    """
    #!/bin/bash
    mkdir -p ${sample_id}
    module load samtools

    samtools index -b ${bam_file} "${sample_id}/${bam_file}.bai"


    """
}


process INDEX_BAM_HIS_TUMOR {
     label "very_short_process"
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
    publishDir "${projectDir}/output/tumor", mode: "copy"

    input:
   tuple val(ix), val(sample_id), path(bam_file)

    output:
    tuple val(ix), val(sample_id), path(bam_file), path("${sample_id}/${bam_file}.bai"), emit: output

    script:
    """
    #!/bin/bash
    mkdir -p ${sample_id}
    module load samtools

    samtools index -b ${bam_file} "${sample_id}/${bam_file}.bai"


    """
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
    tuple val(ix), val(sample_id), path(tumor_hisat_bam), path(tumor_hisat_bam_bai), path(normal_tumor_hisat_bam), path(normal_tumor_hisat_bam_bai), path(snp_mut_intervals), path(snp_mut_intervals_bed)

    output:
    path "${sample_id}/${sample_id}_second_call_stats.*"
    path "${sample_id}/${sample_id}_second_mutect.vcf.*"

    script:
    template "mutect_R2.sh"
}