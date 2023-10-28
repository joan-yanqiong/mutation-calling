# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
cmd_args <- commandArgs(trailingOnly = FALSE)
has_script_filepath <- startsWith(cmd_args, "--file=")
if (sum(has_script_filepath)) {
    setwd(dirname(unlist(strsplit(cmd_args[has_script_filepath], "=")))[2])
}

renv::load(here::here())


# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Formatting mutations from RNAseq",
    )
    parser$add_argument(
        "--input_file",
        type = "character",
        required = TRUE,
        help = "Input file with mutations"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/mutations_postprocessed")
    args$input_file <- glue("{here::here()}/output/tumor/SRR5088829/SRR5088829_second_mutect.annovar.hg19_multianno.txt")
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

annot_mut_files <- list.files(glue("{here::here()}/output/mutations"), pattern = "hg19_multianno.txt", full.names = TRUE)

for (i in annot_mut_files) {
    log_info(glue("Current file: {i}"))
    args$input_file <- i

    log_info("Extract sample id...")
    sample_id <- str_split(get_name(args$input_file), "\\.", simplify = TRUE)[1]
    log_info(glue("Current sample id: {sample_id}"))

    log_info("Load mutations...")
    mutations <- data.frame(fread(args$input_file, sep = "\t")) %>%
        select(-c(V12)) %>%
        rename(allel_freq = V13) %>%
        filter(Func.refGene == "exonic")
    # , ExonicFunc.refGene %in% c("nonsynonymous_SNV", "unknown", "stopgain"))

    log_info("Save file...")
    write.table(mutations, file = glue("{args$output_dir}/{sample_id}_mutations.txt"), sep = "\t", row.names = FALSE)
}
log_info("COMPLETED!")
