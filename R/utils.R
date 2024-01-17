#' Extract name from filepath
#'
#' @param filepath filepath
#' @return name of file without extension
#'
#' @examples examples get_name("R/my_script.R") returns 'my_script'
#' @export
#' @importFrom tools file_path_sans_ext
#' @importFrom fs path_file
#'
get_name <- function(filepath) {
    return(file_path_sans_ext(path_file(filepath)))
}

#' Create directory
#'
#' @description Create directory if it does not exist
#'
#' @param dir_path Path to directory to be created
#'
#' @examples create_dir("/Users/johndoe/my_new_dir")
#' @export
#' @importFrom glue glue
create_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        log_info(glue("Creating directory {dir_path}"))
        dir.create(dir_path, recursive = TRUE)
    } else {
        log_warn(glue("Directory {dir_path} already exists."))
    }
}

#' Formatting VCF files
#' @param file VCF file
#' @return Data frame with formatted VCF file
#' @description Format VCF files for downstream analysis
#' @export
format_mutations <- function(file) {
    curr_vcf <- vcfR::read.vcfR(file, verbose = FALSE, check_keys = TRUE)
    curr_info <- vcfR::extract_info_tidy(curr_vcf)
    curr_gt <- vcfR::extract_gt_tidy(curr_vcf,)
    curr_fix <- vcfR::getFIX(curr_vcf)
    return(merge(curr_gt,
        cbind(curr_info, curr_fix),
        all.x = TRUE
    ))
}