include {TUMOR_REALIGN_PREPROCESS; NORMAL_REALIGN_PREPROCESS; TUMOR_REALIGN_ANNOTATE; NORMAL_REALIGN_ANNOTATE; TUMOR_HISAT_ALIGN; PAIR_HISAT_ALIGN; MUTECT_R2 } from "../nf-modules/rna_mutect2"
include {HISAT_ALIGN_TUMOR_BAM; HISAT_ALIGN_PAIR_BAM} from "../nf-modules/utils"

workflow RNA_MUTECT {
  take:
  mutect_round1
  rna_tumor_processed
  dna_normal_processed
  ref_genome
  index_dir
  dbSNP_set
  cosmic_set

  main:

  pairs_info = Channel.fromPath(params.sample_sheet) \
  | splitCsv(header:true) \
  | map { row -> tuple(row.ix, row.tumor_id, row.normal_id, row.pair_id) }


  // Conversion of BAMs to FASTQs
  // 1a. Create channels for tumor and normal
  realign_input_tumor = pairs_info.join(
    mutect_round1, by: 0
  ).join(rna_tumor_processed, by: 0)

  realign_input_pair = pairs_info.join(
    mutect_round1, by: 0
  ).join(dna_normal_processed, by: 0)

  TUMOR_REALIGN_PREPROCESS(
    jar_files = params.jar_files,
    realign_input_tumor,
    ref_genome
  )
  // 2. Normal
  // tuple(tumor_id, normal_id, pair, maf, maf_idx, recal_bam, recal_bai

  NORMAL_REALIGN_PREPROCESS(
    jar_files = params.jar_files,
    realign_input_pair,
    ref_genome)

  // // ---- REALIGNING WITH HISAT2 ---- //
  // //  ANNOTATE probably redundant
  // // TUMOR_REALIGN_ANNOTATE(TUMOR_REALIGN_PREPROCESS.out.output)
  // // NORMAL_REALIGN_ANNOTATE(NORMAL_REALIGN_PREPROCESS.out.output)

  // tumor_realign_preproc = tuple("SRR5088818", "SRR5134767","SRR5134767_SRR5088818", file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/tmp/SRR5088818_tmp_sequence_1.fastq"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/tmp/SRR5088818_tmp_sequence_2.fastq"), "tumor")
  // TUMOR_HISAT_ALIGN(index_dir = index_dir, tumor_realign_preproc)

  // TUMOR_HISAT_ALIGN(index_dir = index_dir, TUMOR_REALIGN_PREPROCESS.out.output)
  // PAIR_HISAT_ALIGN(index_dir = index_dir, NORMAL_REALIGN_PREPROCESS.out.output)

  // // // Convert SAM to BAM
  // HISAT_ALIGN_TUMOR_BAM(TUMOR_HISAT_ALIGN.out)
  // HISAT_ALIGN_PAIR_BAM(PAIR_HISAT_ALIGN.out)

  // tumor_realign_preprc_mutations = Channel.of(tuple("SRR5088818", file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/tmp/snp_mutations.intervals"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/tmp/snp_mutations.intervals.bed")))


  // // mutect_r2_input = HISAT_ALIGN_TUMOR_BAM.out.join(
  // //   HISAT_ALIGN_PAIR_BAM.out,
  // //   by: 0
  // // ).join(tumor_realign_preprc_mutations, by: 0)


  // mutect_r2_input = HISAT_ALIGN_TUMOR_BAM.out.join(
  //   HISAT_ALIGN_PAIR_BAM.out,
  //   by: 0
  // ).join(TUMOR_REALIGN_PREPROCESS.out.snp_mut_intervals, by: 0)

  // // ---- MUTECT ROUND2 ---- //
  // tuple val(sample_id), path(tumor_hisat_bam), path(normal_tumor_hisat_bam), path(snp_mut_intervals), path(snp_mut_intervals_bed)
  // MUTECT_R2(
  //   ref_genome,
  //   dbSNP_set,
  //   cosmic_set,
  //   mutect_r2_input
  // )
}