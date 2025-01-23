#' @title Cleaning the GWAS summary data for colocalization
#'
#' @description This program is designed to convert GWAS summary data intended for subsequent colocalization analyses into a specific tidy format. Given the GWAS phenotype categories and the information available, users must select one of the following input configurations based on their particular circumstances (see details).
#'
#' @details
#'  In any case, your GWAS summary data should include these following columns: SNP, CHR, BP. In addition, you need to provide extra content based on your study, specifically:
#'
#'  For quantitative GWAS phenotypes, you should:
#'
#'  Specify gwas_data_type = "quant", and then choose one of the following scenarios:
#'
#'    1.With phenotype standard deviation (sdY) and variants' beta and standard error:
#'    Please provide the "gwas_sdy_input" parameter, and ensure the GWAS summary data includes the following columns: beta, SE.
#'    Note: Although the coloc package requires varbeta (i.e., the square of SE), only SE needs to be provided here. The same applies to subsequent cases.
#'
#'    2.With phenotype standard deviation (sdY) but without variants' beta and standard error:
#'    Please provide the "gwas_sdy_input" parameter and ensure the GWAS summary data includes the following columns: pvalues, MAF.
#'
#'    3.Without phenotype standard deviation (sdY) but with variants' beta and standard error:
#'    Please ensure the GWAS summary data includes the following columns: beta, SE, MAF. You should provide the sample size of GWAS in next stage, too.
#'
#'    4.Without phenotype standard deviation (sdY) and without variants' beta and standard error:
#'    Please ensure the GWAS summary data includes the following columns: MAF, pvalues. You should provide the sample size of GWAS in next stage, too.
#'
#'  For binary (Case-Control) GWAS phenotypes, you should:
#'
#'  Specify gwas_data_type = "cc", and provide the "gwas_s_input" parameter (requisite). Then, choose one of the following scenarios:
#'
#'    1.With variants' beta and standard error:
#'    Please ensure nsure the GWAS summary data includes the following columns: beta, SE.
#'
#'    2.Without variants' beta and standard error:
#'    Please ensure the GWAS summary data includes the following columns: pvalues, MAF.
#'
#'  It should be noted that if you wish to align the strands between GWAS and eQTL in the future, you need to set "gwas_allele_input" to TRUE and provide additional column in GWAS summary: A1 (effect allele) and A2 (other allele).
#'
#'
#' @param pheno_name
#' @param SNP_data_create_above
#' @param gwas_data_folder
#' @param gwas_data_suffix
#' @param gwas_allele_input
#' @param gwas_data_type
#' @param gwas_sdy_input
#' @param gwas_beta_input
#' @param save_folder
#'
#' @return
#' @export
#'
#' @examples

gwas_data_preparation <- function(pheno_name, SNP_data_create_above, gwas_data_folder, gwas_data_suffix, gwas_allele_input,
                                  gwas_data_type, gwas_sdy_input = NULL, gwas_s_input = NULL, gwas_beta_input, save_folder)
{
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

  # check SNP_data_create_above
  if (missing(SNP_data_create_above) || !is.character(SNP_data_create_above) || length(SNP_data_create_above) != 1 || !file.exists(SNP_data_create_above)) {
    rlang::abort("'SNP_data_create_above' must be provided as a file path.")
  }
  SNP_data <- data.table::fread(SNP_data_create_above, stringsAsFactors = FALSE, data.table=FALSE)
  if (!is.data.frame(SNP_data) && !is.matrix(SNP_data)) {
    rlang::abort("'SNP_data' must be a valid file path or a data frame/matrix.")
  }
  expected_colnames <- c("SNP", "CHR", "BP", "GENE", "START", "END", "GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "PHENO")
  actual_colnames <- colnames(SNP_data)
  if (!identical(actual_colnames, expected_colnames)) {
    rlang::abort(paste0(
      "The expected SNP_data file column names should be:\n",
      paste(expected_colnames, collapse = "\t"), "\n",
      "But the provided file column names are:\n",
      paste(actual_colnames, collapse = "\t"), "\n",
      "Please check the input file for any issues."
    ))
  }

  if (!all(pheno_name %in% unique(SNP_data$PHENO))) {
    missing_phenos <- pheno_name[!pheno_name %in% unique(SNP_data$PHENO)]
    rlang::abort(paste0("Error: The following phenos are missing in PHENO column of SNP_data: ", paste(missing_phenos, collapse = ", ")))
  }

  ##check gwas_data_suffix
  if (is.null(gwas_data_suffix) || !(is.character(gwas_data_suffix) && length(gwas_data_suffix) == 1)) {
    rlang::abort("Error: 'gwas_data_suffix' must be provided and must be a single string when 'smr_format' is FALSE.")
  }

  ##check gwas_data_folder
  if (missing(gwas_data_folder) || !is.character(gwas_data_folder) || length(gwas_data_folder) != 1 || !dir.exists(gwas_data_folder)) {
    rlang::abort("'gwas_data_folder' must be provided as a file path.")
  }
  gwas_data_folder <- normalizePath(gwas_data_folder)
  for (pheno in pheno_name) {
    file_path <- file.path(gwas_data_folder, paste0(pheno, gwas_data_suffix))
    if (!file.exists(file_path)) {
      rlang::abort(paste0("The GWAS file for pheno '", pheno, "' does not exist: ", file_path))
    }
  }

  # check gwas_allele_input
  if (!is.logical(gwas_allele_input) || length(gwas_allele_input) != 1) {
    rlang::abort("'gwas_allele_input' must be either TRUE or FALSE.")
  }

  # check gwas_data_type
  if (is.null(gwas_data_type) || !(gwas_data_type %in% c("quant", "cc")) || length(gwas_data_type) != 1) {
    rlang::abort("Error: 'gwas_data_type' must be either 'quant' or 'cc'.")
  }

  # check gwas_sdy_input
  if (!is.logical(gwas_sdy_input) && !is.null(gwas_sdy_input)) {
    rlang::abort("Error: 'gwas_sdy_input' must be logical or NULL.")
  }
  if (gwas_data_type == "quant" && is.null(gwas_sdy_input)) {
    rlang::abort("Error: 'gwas_sdy_input' cannot be NULL when 'gwas_data_type' is 'quant'.")
  }
  if (gwas_data_type == "cc" && !is.null(gwas_sdy_input)) {
    rlang::abort("Error: 'gwas_sdy_input' must be NULL when 'gwas_data_type' is 'cc'.")
  }

  # check gwas_s_input
  if (!is.logical(gwas_s_input) && !is.null(gwas_s_input)) {
    rlang::abort("Error: 'gwas_s_input' must be logical or NULL.")
  }
  if (gwas_data_type == "cc" && is.null(gwas_s_input)) {
    rlang::abort("Error: 'gwas_s_input' cannot be NULL when 'gwas_data_type' is 'cc'.")
  }
  if (gwas_data_type == "quant" && !is.null(gwas_s_input)) {
    rlang::abort("Error: 'gwas_s_input' must be NULL when 'gwas_data_type' is 'quant'.")
  }

  # check gwas_beta_input
  if (!is.logical(gwas_beta_input) || length(gwas_beta_input) != 1) {
    rlang::abort("Error: 'gwas_beta_input' must be a logical value.")
  }


  # check save_folder
  if (missing(save_folder) || !is.character(save_folder) || length(save_folder) != 1 || !dir.exists(save_folder)) {
    rlang::abort("'save_folder' must be provided and must be a valid directory path.")
  }
  if (!file.access(save_folder, 2) == 0) {
    rlang::abort("'save_folder' must be a directory where the user has write permissions.")
  }
  save_folder <- normalizePath(save_folder)




  gwas_preparation <- function(pheno, SNP_data, gwas_data_folder, gwas_data_suffix, gwas_allele_input,
                               gwas_data_type, gwas_sdy_input, gwas_beta_input, save_folder)
  {
    gwas_file_path <- file.path(gwas_data_folder, paste0(pheno, gwas_data_suffix))
    gwas_data <- data.table::fread(gwas_file_path, data.table = FALSE, stringsAsFactors = FALSE)
    if (!is.data.frame(gwas_data) && !is.matrix(gwas_data)) {
      rlang::abort("Error: GWAS data must be a data.frame or matrix.")
    }

    colname_map <- list(
      SNP = c("rsid", "variant", "snp", "markername"),
      CHR = c("chrom", "chr"),
      BP = c("position", "basepair", "bp"),
      MAF = c("maf","freq","af1","af","freq1"),
      BETA = c("beta", "b", "effect"),
      SE = c("s.e.", "se", "standard_error", "stderr"),
      P = c("p", "pval", "pvalues", "pvalue"),
      A1 = c("allele1", "effect_allele", "a1"),
      A2 = c("allele0", "allele2", "other_allele")
    )

    # standardized column names
    standardize_colnames <- function(colnames, colname_map) {
      for (standard_col in names(colname_map)) {
        colnames_lower <- tolower(colnames)
        alt_names <- tolower(colname_map[[standard_col]])

        # check colnames
        colnames[colnames_lower %in% alt_names] <- standard_col
      }
      return(colnames)
    }
    colnames(gwas_data) <- standardize_colnames(colnames(gwas_data), colname_map)

    required_cols <- c()

    if (gwas_data_type == "quant") {
      if (gwas_sdy_input == TRUE && gwas_beta_input == TRUE) {
        required_cols <- c("SNP", "CHR", "BP", "BETA", "SE")
      } else if (gwas_sdy_input == FALSE && gwas_beta_input == TRUE) {
        required_cols <- c("SNP", "CHR", "BP", "BETA", "SE", "MAF")
      } else if (gwas_sdy_input == FALSE && gwas_beta_input == FALSE) {
        required_cols <- c("SNP", "CHR", "BP", "P", "MAF")
      } else if (gwas_sdy_input == TRUE && gwas_beta_input == FALSE) {
        required_cols <- c("SNP", "CHR", "BP", "P", "MAF")
      }
    } else if (gwas_data_type == "cc") {
      if (gwas_beta_input == TRUE) {
        required_cols <- c("SNP", "CHR", "BP", "BETA", "SE")
      } else {
        required_cols <- c("SNP", "CHR", "BP", "P", "MAF")
      }
    }
    if (gwas_allele_input == TRUE) {
      required_cols <- c(required_cols, "A1", "A2")
    }

    missing_cols <- setdiff(tolower(required_cols), tolower(colnames(gwas_data)))
    if (length(missing_cols) > 0) {
      rlang::abort(paste0("Error: The following required columns are missing from the GWAS data: ", paste(missing_cols, collapse = ", ")))
    }

    gwas_data <- gwas_data[, required_cols, drop = FALSE]
    if ("MAF" %in% colnames(gwas_data)) {
      gwas_data$MAF[gwas_data$MAF > 0.5] <- 1 - gwas_data$MAF[gwas_data$MAF > 0.5]
    }

    dir.create(file.path(save_folder, pheno), recursive = TRUE)
    for (i in 1:nrow(SNP_data)) {
      snp_row <- SNP_data[i, ]
      chr_value <- snp_row$CHR
      start_value <- snp_row$START
      end_value <- snp_row$END
      filtered_data <- gwas_data[gwas_data$CHR == chr_value & gwas_data$BP >= start_value & gwas_data$BP <= end_value, ]
      output_file <- file.path(save_folder, pheno, paste0(snp_row$SNP, ".txt"))
      write.table(filtered_data, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }
  }

  for (pheno in pheno_name){
    current_SNP_data <- dplyr::filter(SNP_data,PHENO==pheno)
    gwas_preparation(pheno = pheno, SNP_data = current_SNP_data, gwas_data_folder = gwas_data_folder, gwas_data_suffix = gwas_data_suffix, gwas_allele_input = gwas_allele_input,
                                 gwas_data_type = gwas_data_type, gwas_sdy_input = gwas_sdy_input, gwas_beta_input = gwas_beta_input, save_folder = save_folder)
    message(paste0("Data preparation for pheno '", pheno, "' completed successfully!"))
  }
}


