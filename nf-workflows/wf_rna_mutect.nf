include {TUMOR_REALIGN_PREPROCESS; NORMAL_REALIGN_PREPROCESS; TUMOR_HISAT_ALIGN; PAIR_HISAT_ALIGN; MUTECT_R2 ; INDEX_BAM_HIS_NORMAL; INDEX_BAM_HIS_TUMOR; SORT_BAM_COORD_HIS_NORMAL; SORT_BAM_COORD_HIS_TUMOR} from "../nf-modules/rna_mutect2"
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

  // mutect_r1 = Channel.of(tuple("1", "SRR5088818",  file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_onco.maf"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_mutect.vcf.idx")), tuple("2", "SRR5088829",  file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088829/SRR5088829_onco.maf"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088829/SRR5088829_mutect.vcf.idx")))

  // Conversion of BAMs to FASTQs
  // 1a. Create channels for tumor and normal
  realign_input_tumor = pairs_info.join(
    mutect_round1, by: [0,1]
  ).join(rna_tumor_processed, by: [0,1])

  // pairs_info.view()

  // println """
  // pairs-info:${pairs_info}
  // dna_normal_processed:${dna_normal_processed}
  // rna_tumor_processed:${rna_tumor_processed}
  // mutect_r1:${mutect_round1}
  // """
  // dna_normal_processed.view()
  // rna_tumor_processed.view()
  // mutect_round1.view()
  // pairs_info.view()

  realign_input_pair = pairs_info.join(mutect_round1, by: [0, 1]).join(dna_normal_processed, by: 0).map(it -> it.unique())

  realign_input_pair.view()
  // realign_input_pair2 = pairs_info.join(, by: 0)
  // realign_input_tumor.view()
  // realign_input_pair2.view()
  // // Manual
  // realign_input_tumor = tuple(1, "SRR5088818", "SRR5134767", "SRR5134767_SRR5088818", file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_onco.maf"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_mutect.vcf.idx"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_recal.bam"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_recal.bai"))

  // realign_input_pair = Channel.of(tuple(1, "SRR5088818", "SRR5134767", "SRR5134767_SRR5088818", file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_onco.maf"), file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_mutect.vcf.idx"), file("/cluster/projects/gaitigroup/Users/Joan/nf_work/8f/75f7d06a9a74504e6b067f8ddc03f2/SRR5134767/SRR5134767_recal.bam"), file("/cluster/projects/gaitigroup/Users/Joan/nf_work/8f/75f7d06a9a74504e6b067f8ddc03f2/SRR5134767/SRR5134767_recal.bai")))

  TUMOR_REALIGN_PREPROCESS(
    jar_files = params.jar_files,
    realign_input_tumor,
    ref_genome, "tumor"
  )
//   // 2. Normal
//   // tuple(tumor_id, normal_id, pair, maf, maf_idx, recal_bam, recal_bai

  NORMAL_REALIGN_PREPROCESS(
    jar_files = params.jar_files,
    realign_input_pair,
    ref_genome, "normal")

//   // // ---- REALIGNING WITH HISAT2 ---- //
  TUMOR_HISAT_ALIGN(index_dir = index_dir, TUMOR_REALIGN_PREPROCESS.out.output)
  PAIR_HISAT_ALIGN(index_dir = index_dir, NORMAL_REALIGN_PREPROCESS.out.output)

  // Convert SAM to BAM
  HISAT_ALIGN_TUMOR_BAM(TUMOR_HISAT_ALIGN.out)
  HISAT_ALIGN_PAIR_BAM(PAIR_HISAT_ALIGN.out)

  SORT_BAM_COORD_HIS_NORMAL(HISAT_ALIGN_PAIR_BAM.out, "coordinate", "aligned_hisat2_sortedbycoord") | INDEX_BAM_HIS_NORMAL

  SORT_BAM_COORD_HIS_TUMOR(HISAT_ALIGN_TUMOR_BAM.out, "coordinate", "aligned_hisat2_sortedbycoord") | INDEX_BAM_HIS_TUMOR

  // mutect_r2_input = HISAT_ALIGN_TUMOR_BAM.out.join(
  //   HISAT_ALIGN_PAIR_BAM.out,
  //   by: 0
  //   ).map(it ->
  //   it.unique()).join(TUMOR_REALIGN_PREPROCESS.out.snp_mut_intervals, by:
  //   0).map(it -> it.unique())




  mutect_r2_input = pairs_info.join(INDEX_BAM_HIS_TUMOR.out.output, by: [0, 1]).join(INDEX_BAM_HIS_NORMAL.out.output, by: 0).map(it -> it.unique()).join(TUMOR_REALIGN_PREPROCESS.out.snp_mut_intervals, by: [0, 1])


//   mutect_r2_input = Channel.of(tuple(
//     1,
//     "SRR5088818",
//     file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_hisat_coord.bam"),
//       file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5088818_hisat_coord.bam.bai"),
//     file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5134767_SRR5088818_hisat_coord.bam"),
//       file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/SRR5134767_SRR5088818_hisat_coord.bam.bai"),
//     file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/snp_mutations.intervals"),
//     file("/cluster/projects/gaitigroup/Users/Joan/h4h-mutation-calling/output/tumor/SRR5088818/snp_mutations.intervals.bed")
//   ))

  // ---- MUTECT ROUND2 ---- //
  MUTECT_R2(
    ref_genome,
    dbSNP_set,
    cosmic_set,
mutect_r2_input
  )
}