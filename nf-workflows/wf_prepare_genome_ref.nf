include { CREATE_FASTA_DICT; CREATE_FASTA_INDEX } from "../nf-modules/prep-ref-genome"

workflow SET_REF_GENOME {
    main:
    ref_genome_basepath = "${file(params.ref_genome).parent}/${file(params.ref_genome).simpleName}"
    fai_path = file("${ref_genome_basepath}.fasta.fai")
    dict_path = file("${ref_genome_basepath}.dict")

    fasta_exists = file("${params.ref_genome}").exists()
    fai_exists = fai_path.exists()
    dict_exists = dict_path.exists()

    if (fasta_exists && fai_exists && dict_exists) {
        println "FASTA, FAI, and DICT files already exist"
        ref_genome = tuple(params.ref_genome, file("${ref_genome_basepath}.fasta.fai"), file("${ref_genome_basepath}.dict"))
    } else {
        CREATE_FASTA_DICT(params.ref_genome)
        CREATE_FASTA_INDEX(params.ref_genome)
        ref_genome = CREATE_FASTA_DICT.out.join(CREATE_FASTA_INDEX.out, by: 0)
    }
    println ref_genome
    emit:
    ref_genome = ref_genome
}