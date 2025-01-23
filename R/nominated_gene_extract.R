#' Extract genes for colocalization-based linking (i.e., nominated genes)
#'
#' Use this function to organize the list of potentially relevant genes obtained from locus-based linking (cS2G) and similarity-based linking (PoPs)—referred to as “nominated genes”—for subsequent colocalization-based linking work.
#'
#' @param pheno_name
#' @param cs2g_result_folder
#' @param pops_result_folder
#' @param save_folder
#'
#' @return
#' @export
#'
#' @examples

nominated_gene_extract <- function(pheno_name, cs2g_result_folder, pops_result_folder, save_folder)
{
  suppressMessages(suppressWarnings(library(dplyr)))
  # check pheno_name
  check_and_standardize_pheno_name <- function(pheno_name) {
    if (missing(pheno_name) || !is.character(pheno_name) || length(pheno_name) != 1 || nchar(pheno_name) == 0) {
      rlang::abort("'pheno_name' must be provided and must be a non-empty character string.")
    }
    pheno_name <- unlist(strsplit(pheno_name, ","))
    pheno_name <- trimws(pheno_name)
    if (any(nchar(pheno_name) == 0)) {
      rlang::abort("'pheno_name' contains empty values after splitting by commas.")
    }
    return(pheno_name)
  }
  pheno_name <- check_and_standardize_pheno_name(pheno_name)

  if (missing(cs2g_result_folder) || !is.character(cs2g_result_folder) || length(cs2g_result_folder) != 1 || !dir.exists(cs2g_result_folder)) {
    rlang::abort("'cs2g_result_folder' must be provided and must be a valid directory path.")
  }
  cs2g_result_folder <- normalizePath(cs2g_result_folder)

  if (missing(pops_result_folder) || !is.character(pops_result_folder) || length(pops_result_folder) != 1 || !dir.exists(pops_result_folder)) {
    rlang::abort("'pops_result_folder' must be provided and must be a valid directory path.")
  }
  pops_result_folder <- normalizePath(pops_result_folder)

  if (missing(save_folder) || !is.character(save_folder) || length(save_folder) != 1 || !dir.exists(save_folder)) {
    rlang::abort("'save_folder' must be provided and must be a valid directory path.")
  }
  if (!file.access(save_folder, 2) == 0) {
    rlang::abort("'save_folder' must be a directory where the user has write permissions.")
  }
  save_folder <- normalizePath(save_folder)

  gene_extract <- function(pheno, cs2g_result_folder, pops_result_folder, save_folder)
  {
    cs2g_data <- data.table::fread(file.path(cs2g_result_folder,paste0(pheno,"_cS2G.csv")), data.table=F, stringsAsFactors = F)
    if (!is.data.frame(cs2g_data) && !is.matrix(cs2g_data)) {
      rlang::abort("'cs2g_data' must be a valid file path or a data frame/matrix.")
    }
    expected_colnames <- c("SNP", "GENE", "VALUE", "STRATEGY", "VALUE_WEIGHTED", "precision", "PIP", "confidence_score", "ENSGID")
    actual_colnames <- colnames(cs2g_data)
    if (!identical(actual_colnames, expected_colnames)) {
      rlang::abort(paste0(
        "The expected pops_result file column names should be:\n",
        paste(expected_colnames, collapse = "\t"), "\n",
        "But the provided file column names are:\n",
        paste(actual_colnames, collapse = "\t"), "\n",
        "Please check the input file for any issues."
      ))
    }

    pops_data <- data.table::fread(file.path(pops_result_folder,paste0(pheno,"_pops_result_processed.csv")), data.table=F, stringsAsFactors = F)
    if (!is.data.frame(pops_data) && !is.matrix(pops_data)) {
      rlang::abort("'pops_data' must be a valid file path or a data frame/matrix.")
    }
    expected_colnames <- c("CHR", "SNP", "BP", "GENE", "START", "END", "GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
    actual_colnames <- colnames(pops_data)
    if (!identical(actual_colnames, expected_colnames)) {
      rlang::abort(paste0(
        "The expected pops_result file column names should be:\n",
        paste(expected_colnames, collapse = "\t"), "\n",
        "But the provided file column names are:\n",
        paste(actual_colnames, collapse = "\t"), "\n",
        "Please check the input file for any issues."
      ))
    }

    gene_list <- unique(na.omit(c(cs2g_data$ENSGID[cs2g_data$confidence_score < 0.8],
                                  unlist(pops_data[, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")])) %>% .[. != "NONE"]))

    merged_data <- merge(pops_data, cs2g_data[cs2g_data$confidence_score < 0.8, c("SNP", "ENSGID")], by = "SNP", all.x = TRUE)
    merged_data$GENE6 <- merged_data$ENSGID
    merged_data$GENE6[is.na(merged_data$GENE6)] <- "NONE"
    merged_data$ENSGID <- NULL
    merged_data$PHENO <- pheno

    cs2g_data_cred <- subset(cs2g_data,cs2g_data$confidence_score >= 0.8)
    cs2g_data_cred$PHENO <- pheno

    return(list(gene_list = gene_list, merged_data = merged_data, cs2g_data_cred = cs2g_data_cred))
  }

  combined_gene_list <- c()
  combined_merged_data <- data.frame()
  combined_cs2g_data_cred <- data.frame()

  for (pheno in pheno_name) {
    result_list <- gene_extract(pheno, cs2g_result_folder, pops_result_folder, save_folder)
    combined_gene_list <- c(combined_gene_list, result_list$gene_list)
    combined_merged_data <- rbind(combined_merged_data, result_list$merged_data)
    combined_cs2g_data_cred <- rbind(combined_cs2g_data_cred, result_list$cs2g_data_cred)
  }
  write.table(unique(combined_gene_list), file = file.path(save_folder, "gene_for_coloc.txt"),row.names = FALSE, col.names = FALSE, sep="\t", quote=F)
  write.csv(combined_merged_data, file = file.path(save_folder, "SNP_data_for_coloc.csv"), row.names=FALSE)
  write.csv(combined_cs2g_data_cred, file = file.path(save_folder, "cs2g_cred_gene.csv"), row.names=FALSE)
}
