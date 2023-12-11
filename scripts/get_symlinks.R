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
    args$output_dir <- glue("{here::here()}/output/")
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
file_paths <- list.files("/Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-mutation-calling/output/test_set", recursive = TRUE, full.names = TRUE)

true_file_paths <- sapply(file_paths, Sys.readlink)
df <- data.frame(true_file_paths) %>%
    rownames_to_column("file_paths") %>%
    filter(true_file_paths != "")
# all_file_paths <- cbind(file_paths, true_file_paths)
# true_file_paths <- true_file_paths[true_file_paths != ""]

write.table(df, file = "output/tumor_true_file_paths.txt", row.names = FALSE, quote = FALSE, sep = "\t")
