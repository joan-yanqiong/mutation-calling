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
    args$output_dir <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/000_misc/Lupus"
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
log_info("Get list of downloaded samples...")
log_info("Check RNAseq samples...")
rnaseq_condition <- list.files("/Users/joankant/Desktop/gaitigroup/Users/Joan/001_data/Lupus/Riaz/rnaseq_condition")
log_info(glue("Number of RNAseq samples: {length(rnaseq_condition)}"))

log_info("Check WES samples...")
wes_condition <- list.files("/Users/joankant/Desktop/gaitigroup/Users/Joan/001_data/Lupus/Riaz/wes_condition")
log_info(glue("Number of WES samples: {length(wes_condition)}"))

wes_normal <- list.files("/Users/joankant/Desktop/gaitigroup/Users/Joan/001_data/Lupus/Riaz/wes_normal")
log_info(glue("Number of WES normal samples: {length(wes_normal)}"))

rnaseq_condition_df <- cbind(rnaseq_condition, c(rep("RNAseq_condition", length(rnaseq_condition))))

wes_condition_df <- cbind(wes_condition, c(rep("WES_condition", length(wes_condition))))
wes_normal_df <- cbind(wes_normal, c(rep("WES_normal", length(wes_normal))))

log_info("Create unified list of downloaded samples...")
downloaded_files <- do.call(rbind, list(wes_condition_df, rnaseq_condition_df, wes_normal_df)) %>%
    data.frame() %>%
    rename(sample_id = wes_condition, group = V2)

date_of_checking <- format(Sys.time(), "%Y%m%d")
write.csv(downloaded_files, glue("{args$output_dir}/riaz_downloaded_samples_{date_of_checking}.csv"), row.names = FALSE)
log_info("COMPLETED!")
