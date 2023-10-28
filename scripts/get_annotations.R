# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

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
        description = "Get metadata",
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/data/sites_to_remove_dbs")
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
log_info("Loading additional libraries...")
pacman::p_load("biomaRt")

log_info("Set up SSL config...")
httr::set_config(httr::config(ssl_verifypeer = FALSE))

log_info("Getting annotations from Ensembl...")
grch37 <- useEnsembl(biomart = "ensembl", GRCh = 37, dataset = "hsapiens_gene_ensembl")
reference <- getBM(attributes = c(
    "ensembl_gene_id",
    "ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "transcript_biotype"
), mart = grch37)

log_info("Keeping only chromosomes of interest (1-22, X, Y)")
tmp <- reference %>% filter(!startsWith(chromosome_name, "HG"), !startsWith(chromosome_name, "CHR"), !startsWith(chromosome_name, "GL"), !startsWith(chromosome_name, "JH"), !startsWith(chromosome_name, "MT"), !startsWith(chromosome_name, "HSCHR"))

log_info("Filtering for pseudo genes")
pseudogenes_2 <- tmp %>%
    filter(str_detect(transcript_biotype, "pseudogene")) %>%
    dplyr::select(chromosome_name, start_position, end_position)
colnames(pseudogenes_2) <- c("chrom", "chromStart", "chromEnd")

log_info("Writing output files...")
write.table(pseudogenes_2, file = glue("{args$output_dir}/pseudo_genes_2.bed"), sep = "\t", row.names = FALSE, quote = FALSE)

log_info("COMPLETED!")