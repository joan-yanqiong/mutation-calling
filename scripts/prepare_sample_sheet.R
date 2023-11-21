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
        description = "Prepare Samplesheet for pipeline", default_output = "misc"
    )
    parser$add_argument(
        "--gseo",
        type = "character",
        required = TRUE,
        help = "GEO metadata file"
    )
    parser$add_argument(
        "--wes_sra",
        type = "character",
        required = TRUE,
        help = "SRA metadata file for WES"
    )
    parser$add_argument(
        "--rnaseq_sra",
        type = "character",
        required = TRUE,
        help = "SRA metadata file for RNAseq"
    )
    parser$add_argument(
        "--fastq_dir",
        type = "character",
        required = TRUE,
        help = "Directory with fastq files"
    )
    parser$add_argument("--normal_samples_oi",
        type = "character",
        required = FALSE,
        help = "File with normal samples of interest", default = NULL
    )
    parser$add_argument("--tumor_samples_oi",
        type = "character",
        required = FALSE,
        help = "File with tumor samples of interest", default = NULL
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("~/Desktop/gaitigroup/Users/Joan/000_misc/Lupus")
    args$gseo <- glue("~/Desktop/gaitigroup/Users/Joan/000_misc/Lupus/GSE91061_series_matrix.xlsx")
    args$fastq_dir <- "/cluster/projects/gaitigroup/Users/Joan/001_data/Lupus/Riaz"
    args$wes_sra <- glue("~/Desktop/gaitigroup/Users/Joan/000_misc/Lupus/SraRunTable_WES.txt")
    args$rnaseq_sra <- glue("~/Desktop/gaitigroup/Users/Joan/000_misc/Lupus/SraRunTable_RNA.txt")
    # args$normal_samples_oi <- glue("/Users/joankant/Library/CloudStorage/OneDrive-UHN/004_Projects/Lupus/pilot_normal_samples_for_download.txt")
    # args$tumor_samples_oi <- NULL
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
log_info("Load additional libraries...")
pacman::p_load(readxl)

log_info("Read metadata from GEO...")
rnaseq_meta <- read_excel(args$gseo, skip = 35)

# Lookup table to match RNAseq with WES
# Columns: rnasq_patient_id, GSM_id
log_info("Create lookup table...")
lookup_table <- data.frame(t(rnaseq_meta[1, ])) %>% rownames_to_column("patient_id")
colnames(lookup_table) <- c("rna_sample_id", "GSM_id")
lookup_table <- lookup_table[-1, ]
head(lookup_table)

# Downloaded from SRA selectors
log_info("Read metadata from SRA selectors...")
wes_meta_sra_normal <- read.csv(args$wes_sra) %>%
    select(Sample.Name, Run) %>%
    rename(wes_normal_id = Run, norm_id = Sample.Name) %>%
    mutate(norm_id = tolower(norm_id)) %>%
    filter(str_detect(norm_id, "norm"))
wes_meta_sra_condition <- read.csv(args$wes_sra) %>%
    select(Sample.Name, Run) %>%
    rename(wes_condition_id = Run, pt_condition_id = Sample.Name) %>%
    mutate(pt_condition_id = tolower(pt_condition_id)) %>%
    filter(!str_detect(pt_condition_id, "norm"))

rnaseq_meta_sra <- read.csv(args$rnaseq_sra) %>%
    select(Sample.Name, Run) %>%
    rename(GSM_id = Sample.Name, rnaseq_condition_id = Run)


# Combine RNAseq and lookup table
log_info("Match metadata RNAseq and lookup table...")
rnaseq_meta_sra <- rnaseq_meta_sra %>%
    left_join(lookup_table, by = "GSM_id") %>%
    mutate(patient_id = str_split(rna_sample_id, simplify = TRUE, pattern = "_")[, 1]) %>%
    mutate(norm_id = tolower(paste0(patient_id, "_norm")), condition = tolower(str_split(rna_sample_id, "_", simplify = TRUE)[, 2])) %>%
    mutate(pt_condition_id = tolower(paste0(patient_id, "_", condition)))

# Combine RNAseq and WES
log_info("Match metadata RNAseq and WES...")
sample_sheet <- rnaseq_meta_sra %>%
    left_join(wes_meta_sra_normal, by = c(norm_id = "norm_id")) %>%
    left_join(wes_meta_sra_condition, by = "pt_condition_id") %>%
    filter(wes_normal_id != "<NA>", wes_condition_id != "<NA>", rnaseq_condition_id != "<NA>")

log_info("Adding sample paths...")
sample_sheet <- sample_sheet %>%
    mutate(
        ix = row_number(),
        normal_tumor_id = paste0(wes_normal_id, "_", rnaseq_condition_id),
        wes_pair_id = paste0(wes_normal_id, "_", wes_condition_id),
        tumor_fastq = paste0(args$fastq_dir, "/rnaseq_condition/", rnaseq_condition_id),
        normal_fastq = paste0(args$fastq_dir, "/wes_normal/", wes_normal_id),
        wes_condition_fastq = paste0(args$fastq_dir, "/wes_condition/", wes_condition_id)
    )

log_info("Save sample sheet...")
write.table(sample_sheet, file = glue("{args$output_dir}/sample_sheet.csv"), sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)
 log_info("COMPLETED!")
