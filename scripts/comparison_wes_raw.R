# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# renv::load(glue(
#     "/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-mutation-calling/renv"
# ))
# Set working directory
cmd_args <- commandArgs(trailingOnly = FALSE)
has_script_filepath <- startsWith(cmd_args, "--file=")
if (sum(has_script_filepath)) {
    setwd(dirname(unlist(strsplit(cmd_args[has_script_filepath], "=")))[2])
}

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Compare to WES",
    )
    parser$add_argument("--wes_mutations", type = "character", required = TRUE, help = "Input file with mutations from WES")
    parser$add_argument("--rna_mutations", type = "character", required = TRUE, help = "Input file with mutations from RNAseq")
    parser$add_argument("--sample_sheet", type = "character", required = TRUE, help = "Sample sheet with patient IDs")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$wes_mutations <- "/Users/joankant/Library/CloudStorage/OneDrive-UHN/004_Projects/Lupus/NIHMS907788-supplement-10.xlsx"
    # args$rna_mutations <- glue("{here::here()}/output/mutations/SRR5088818_mutations.txt")
    args$output_dir <- glue("{here::here()}/output/figures")
    args$sample_sheet <- glue("{here::here()}/misc/sample_sheet.csv")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
log_info("Load libraries...")
pacman::p_load(readxl, VennDiagram)
pacman::p_load_gh("jmw86069/colorjam")

mutation_files <- list.files(glue("{here::here()}/output/mutations_postprocessed"), full.names = TRUE)

log_info("Load WES mutations (reference)...")
wes_mutations <- read_excel(args$wes_mutations, skip = 2) %>% mutate(found_in_dna = 1) #%>% filter(!str_detect(`Variant Classification`, "intron"))
for (i in mutation_files) {
    log_info(glue("Current file: {i}"))
    args$rna_mutations <- i

    log_info("Load RNAseq mutations...")
    rna_mutations <- fread(args$rna_mutations, header = TRUE)
    rna_mutations <- data.frame(rna_mutations) %>% mutate(found_in_rna = 1)

    log_info("Load sample sheet...")
    sample_sheet <- read.table(args$sample_sheet, sep = ",", header = TRUE)

    current_tumor_id <- str_split(basename(args$rna_mutations), "\\.", simplify = TRUE)[1]

    current_patient_id <- sample_sheet %>%
        filter(tumor_id == current_tumor_id) %>%
        pull(patient_id)
    log_info(glue("Current patient ID: {current_patient_id} with tumor ID: {current_tumor_id}"))

    log_info("Filtering mutations...")
    patient_mutations <- wes_mutations %>% filter(Patient == current_patient_id)

    rna_mutations_prepped <- rna_mutations %>%
        rename(Chromosome = Chr) %>%
        mutate(mutation = paste0(Ref, ">", Alt)) %>%
        select(Chromosome, mutation, found_in_rna, allel_freq, Start, End, Func.refGene, Gene.refGene)
    patient_mutations_prepped <- patient_mutations %>%
        mutate(mutation = str_sub(HGVS_c, start = -3)) %>%
        select(Chromosome, mutation, `Hugo Symbol`, found_in_dna, Start, End, `Variant Classification`)

    merge_ids <- intersect(colnames(rna_mutations_prepped), colnames(patient_mutations_prepped))

    log_info(glue("Number of mutations in RNAseq: {nrow(rna_mutations_prepped)}"))
    log_info(glue("Number of mutations in WES: {nrow(patient_mutations_prepped)}"))
    combined <- merge(rna_mutations_prepped, patient_mutations_prepped, by = merge_ids, all = TRUE)

    log_info(glue("Total number of mutations: {nrow(combined)}"))

    combined <- combined %>%
        mutate(common = found_in_dna + found_in_rna) %>%
        mutate(mut_id = paste0("Chr_", Chromosome, "_", Start, "_", End, ":", mutation))
    common_interactions <- combined %>% filter(common == 2)

    log_info("Save mutations...")
    write.table(combined, file = glue(
        "{args$output_dir}/{current_patient_id}_compare_mutations.csv"
    ))

    log_info(glue("Number of common mutations: {nrow(common_interactions)} ({nrow(common_interactions) / nrow(patient_mutations_prepped) * 100}%)"))

    log_info("Create VennDiagrams...")
    bulk_rnaseq <- combined %>% pull(found_in_rna)
    bulk_rnaseq[is.na(bulk_rnaseq)] <- 0
    wes <- combined %>% pull(found_in_dna)
    wes[is.na(wes)] <- 0
    x <- list(bulk_rnaseq = combined %>% filter(found_in_rna == 1) %>% pull(mut_id), wes = combined %>% filter(found_in_dna == 1) %>% pull(mut_id))
    category.names <- names(x)
    filename <- glue("{args$output_dir}/{current_tumor_id}_venn.png")

    create_venndiagram(x, category.names, filename, main = glue("Tumor ID: {current_tumor_id}"))
}
log_info("COMPLETED!")
