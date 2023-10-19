include { BWA_INDEX; STAR_INDEX; HISAT_INDEX } from "../nf-modules/build_indices"

workflow BUILD_INDICES {
    take:
        ref_genome

    main:
    if (file(params.bwa_index_dir).isEmpty()) {
        BWA_INDEX(
        index_dir = params.bwa_index_dir,
        ref_path = ref_genome
    )
        bwa_index = BWA_INDEX.out.indexed_dir
    } else {
        println "BWA index already exists, skipping"
        bwa_index = file(params.bwa_index_dir)
    }

    if (file(params.star_index_dir).isEmpty()) {
        STAR_INDEX(
        index_dir = params.star_index_dir,
        ref_path = ref_genome,
        params.gtf_path
        )
        star_index = STAR_INDEX.out.indexed_dir
    } else {
        println "STAR index already exists, skipping"
        star_index = file(params.star_index_dir)
    }

    if (file(params.hisat_index_dir).isEmpty()) {
        HISAT_INDEX(
            index_dir = params.hisat_index_dir,
            ref_path = ref_genome
        )
        hisat_index = HISAT_INDEX.out.indexed_dir
    } else {
        println "HISAT index already exists, skipping"
        hisat_index = file(params.hisat_index_dir)
    }

    println """\
    =======================
    ---- BUILD INDICES ----
    =======================
    BWA index: ${bwa_index}
    STAR index: ${star_index}
    HISAT index: ${hisat_index}
    """.stripIndent()

    emit:
    hisat_index = hisat_index
    star_index = star_index
    bwa_index = bwa_index
}