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
    args$output_dir <- "~/Desktop/gaitigroup/Users/Joan/000_misc/Lupus"
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
# sample_sheet <- read.csv("/Users/joankant/Desktop/gaitigroup/Users/Joan/000_misc/Lupus/sample_sheet.csv")


# sample_sheet_sub <- read.csv("/Users/joankant/Desktop/gaitigroup/Users/Joan/000_misc/Lupus/riaz_wes_test_set.csv")

# wes_cond_id <- sample_sheet_sub %>% pull(wes_cond_id)

# sample_sheet_sub_updated <- sample_sheet %>%
#     filter(wes_condition_id %in% wes_cond_id) %>%
#     mutate(ix = row_number())

out <- read.csv("/Users/joankant/Desktop/gaitigroup/Users/Joan/000_misc/Lupus/test_set.csv")

write.table(out, file = glue("{args$output_dir}/test_set.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
