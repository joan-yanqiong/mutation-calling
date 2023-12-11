include {TUMOR_REALIGN_PREPROCESS; NORMAL_REALIGN_PREPROCESS; TUMOR_HISAT_ALIGN; PAIR_HISAT_ALIGN; MUTECT_R2 ; INDEX_BAM_HIS_NORMAL; INDEX_BAM_HIS_TUMOR; SORT_BAM_COORD_HIS_NORMAL; SORT_BAM_COORD_HIS_TUMOR;FILTERING; FUNCOTATOR; ANNOVAR_R2} from "../nf-modules/rna_mutect2"
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
  | map { row -> tuple(row.ix, row.rnaseq_tumor_id, row.wes_normal_id, row.normal_tumor_id) }
  normal_pairs_info = Channel.fromPath(params.sample_sheet) \
  | splitCsv(header:true) \
  | map { row -> tuple(row.wes_normal_id, row.ix, row.wes_normal_id) }

  pairs_info.view()
  // normal_pairs_info.view()

  dna_normal_processed_adapted = normal_pairs_info.combine(dna_normal_processed, by:0).map( it -> it.subList(1, it.size()))

  dna_normal_processed_adapted.view()

  // Conversion of BAMs to FASTQs
  // 1a. Create channels for tumor and normal
  realign_input_tumor = pairs_info.join(
    mutect_round1, by: [0,1]
  ).join(rna_tumor_processed, by: [0,1])


  realign_input_pair = pairs_info.join(mutect_round1, by: [0, 1]).join(dna_normal_processed_adapted, by: 0).map(it -> it.unique())

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
  HISAT_ALIGN_TUMOR_BAM(TUMOR_HISAT_ALIGN.out.output)
  HISAT_ALIGN_PAIR_BAM(PAIR_HISAT_ALIGN.out.output)

  SORT_BAM_COORD_HIS_NORMAL(HISAT_ALIGN_PAIR_BAM.out, "coordinate")

  INDEX_BAM_HIS_NORMAL(SORT_BAM_COORD_HIS_NORMAL.out.output)

  SORT_BAM_COORD_HIS_TUMOR(HISAT_ALIGN_TUMOR_BAM.out, "coordinate")

  INDEX_BAM_HIS_TUMOR(SORT_BAM_COORD_HIS_TUMOR.out.output)

  mutect_r2_input = pairs_info.join(INDEX_BAM_HIS_TUMOR.out.output, by: [0, 1]).join(INDEX_BAM_HIS_NORMAL.out.output, by: 0).map(it -> it.unique()).join(TUMOR_REALIGN_PREPROCESS.out.snp_mut_intervals, by: [0, 1])

  // ---- MUTECT ROUND2 ---- //
  MUTECT_R2(
    ref_genome,
    dbSNP_set,
    cosmic_set,
    mutect_r2_input,
    suffix = "mutect_R2"
  )

  FILTERING(
    MUTECT_R2.out.output,
    exac_bed = params.exac_bed,
    darned_bed = params.darned_bed,
    radar_bed = params.radar_bed,
    pseudo_genes_bed = params.pseudo_genes_bed,
    pseudo_genes_bed2 = params.pseudo_genes_bed2,
    min_alt_counts = params.min_alt_counts)

  ANNOVAR_R2(
    FILTERING.out.output,
    human_db = params.human_db, suffix = "annovar_R2")

  FUNCOTATOR(
    FILTERING.out.output,
    data_source_path = params.func_data_path,
    suffix = "R2",
    ref_genome = ref_genome,
  )
}
