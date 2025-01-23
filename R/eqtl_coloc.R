#' @title Colocalization-based linking
#'
#' @description This function utilizes "coloc::coloc.abf" to perform colocalization analysis between eQTL and GWAS, providing colocalization-based linking evidence.
#'
#' @param pheno_name
#' @param SNP_data_create_above
#' @param gwas_separate_data_folder
#' @param eqtl_separate_data_folder
#' @param strand_harmonize
#' @param gwas_data_type
#' @param gwas_sdy_input
#' @param gwas_sdy_data
#' @param gwas_s_data
#' @param gwas_beta_input
#' @param gwas_n_data
#' @param eqtl_n_data
#' @param save_folder
#' @param work_part
#' @param coloc_p1
#' @param coloc_p2
#' @param coloc_p12
#'
#' @return
#' @export
#'
#' @examples

eqtl_coloc <- function(pheno_name, SNP_data_create_above, gwas_separate_data_folder, eqtl_separate_data_folder, strand_harmonize, gwas_data_type,
                       gwas_sdy_input, gwas_sdy_data = NULL, gwas_s_data = NULL, gwas_beta_input, gwas_n_data = NULL, eqtl_n_data = NULL, save_folder, work_part = NULL,
                       coloc_p1 = NULL, coloc_p2 = NULL, coloc_p12 = NULL)
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


  # check gwas separate file
  for (pheno in pheno_name) {
    folder_path <- file.path(gwas_separate_data_folder, pheno)

    if (!dir.exists(folder_path)) {
      rlang::abort(paste0("Error: The gwas_separate_data_folder for pheno '", pheno, "' does not exist at path: ", folder_path))
    }
  }
  gwas_separate_data_folder <- normalizePath(gwas_separate_data_folder)

  if (!dir.exists(eqtl_separate_data_folder)) {
    rlang::abort(paste0("Error: The eqtl_separate_data_folder '", eqtl_separate_data_folder, "' does not exist."))
  }
  num_subfolders <- length(list.dirs(eqtl_separate_data_folder, recursive=F,full.names=F))
  message(paste0("There are total ", num_subfolders,
                 " subfolders in '", eqtl_separate_data_folder,
                 "' (always refer to the number of tissues in eqtl), please check whether it matches your expectation."))
  eqtl_separate_data_folder <- normalizePath(eqtl_separate_data_folder)
  tissue_list <- list.dirs(eqtl_separate_data_folder, recursive=F,full.names=F)

  ##check strand_harmonize
  if (!is.logical(strand_harmonize) || length(strand_harmonize) != 1) {
    rlang::abort("Error: 'strand_harmonize' must be a logical value.")
  }

  ##check gwas_data_type
  if (is.null(gwas_data_type) || !(gwas_data_type %in% c("quant", "cc")) || length(gwas_data_type) != 1) {
    rlang::abort("Error: 'gwas_data_type' must be either 'quant' or 'cc'.")
  }

  ##check gwas_sdy_input
  if (!is.logical(gwas_sdy_input) && !is.null(gwas_sdy_input)) {
    rlang::abort("Error: 'gwas_sdy_input' must be logical or NULL.")
  }
  if (gwas_data_type == "quant" && is.null(gwas_sdy_input)) {
    rlang::abort("Error: 'gwas_sdy_input' cannot be NULL when 'gwas_data_type' is 'quant'.")
  }
  if (gwas_data_type == "cc" && !is.null(gwas_sdy_input)) {
    rlang::abort("Error: 'gwas_sdy_input' must be NULL when 'gwas_data_type' is 'cc'.")
  }

  ##check gwas_sdy_data
  if (isTRUE(gwas_sdy_input)) {
    if (is.null(gwas_sdy_data) || !is.character(gwas_sdy_data) || length(gwas_sdy_data) != 1 || !file.exists(gwas_sdy_data)) {
      rlang::abort("Error: 'gwas_sdy_data' must be a non-null character pointing to the path of sdY data file when 'gwas_sdy_input' is TRUE.")
    }
    gwas_sdy_data <- data.table::fread(gwas_sdy_data, data.table=F, stringsAsFactors = F)
    if (!is.data.frame(gwas_sdy_data) && !is.matrix(gwas_sdy_data)) {
      rlang::abort("'gwas_sdy_data' must be a valid file path or a data frame/matrix.")
    }
    colnames_upper <- toupper(colnames(gwas_sdy_data))
    if (!all(c("PHENO", "SDY") %in% colnames_upper)) {
      rlang::abort("Error: gwas_sdy_data must contain columns 'PHENO' and 'SDY' (case-insensitive).")
    }
    colnames(gwas_sdy_data) <- colnames_upper
    gwas_sdy_data <- gwas_sdy_data[, c("PHENO", "SDY")]
    gwas_sdy_data <- setNames(gwas_sdy_data$SDY, gwas_sdy_data$PHENO)
  } else if (is.null(gwas_sdy_input) || gwas_sdy_input == FALSE) {
    if (!is.null(gwas_sdy_data)) {
      rlang::abort("Error: 'gwas_sdy_data' must be NULL when 'gwas_sdy_input' is NULL or FALSE.")
    }
  }

  ##check gwas_s_data
  if (gwas_data_type == "cc") {
    if (is.null(gwas_s_data) || !is.character(gwas_s_data) || length(gwas_s_data) != 1 || !file.exists(gwas_s_data)) {
      rlang::abort("Error: 'gwas_sdy_data' must be a non-null character pointing to the path of s (proportion of case in case-control GWA study) data file when 'gwas_data_type' is cc.")
    }
    gwas_s_data <- data.table::fread(gwas_s_data, data.table=F, stringsAsFactors = F)
    if (!is.data.frame(gwas_s_data) && !is.matrix(gwas_s_data)) {
      rlang::abort("'gwas_s_data' must be a valid file path or a data frame/matrix.")
    }
    colnames_upper <- toupper(colnames(gwas_s_data))
    if (!all(c("PHENO", "S") %in% colnames_upper)) {
      rlang::abort("Error: gwas_s_data must contain columns 'PHENO' and 'S' (case-insensitive).")
    }
    colnames(gwas_s_data) <- colnames_upper
    gwas_s_data <- gwas_s_data[, c("PHENO", "S")]
    gwas_s_data <- setNames(gwas_s_data$S, gwas_s_data$PHENO)
  } else {
    if (!is.null(gwas_s_data)) {
      rlang::abort("Error: 'gwas_s_data' must be NULL when 'gwas_data_type' is quant.")
    }
  }

  ##check gwas_beta_input
  if (!is.logical(gwas_beta_input) || length(gwas_beta_input) != 1) {
    rlang::abort("Error: 'gwas_beta_input' must be a logical value.")
  }

  ##check gwas_n_data
  if (gwas_data_type == "quant" && gwas_sdy_input == FALSE) {
    if (is.null(gwas_n_data) || !is.character(gwas_n_data) || length(gwas_n_data) != 1 || !file.exists(gwas_n_data)) {
      rlang::abort("Error: 'gwas_n_data' must be a non-null character string of length 1 when 'gwas_data_type' is 'quant' and 'gwas_sdy_input' is FALSE.")
    }
    gwas_n_data <- data.table::fread(gwas_n_data, data.table=F, stringsAsFactors = F)
    if (!is.data.frame(gwas_n_data) && !is.matrix(gwas_n_data)) {
      rlang::abort("'gwas_n_data' must be a valid file path or a data frame/matrix.")
    }
    colnames_upper <- toupper(colnames(gwas_n_data))
    if (!all(c("PHENO", "N") %in% colnames_upper)) {
      rlang::abort("Error: gwas_n_data must contain columns 'PHENO' and 'N' (case-insensitive).")
    }
    colnames(gwas_n_data) <- colnames_upper
    gwas_n_data <- gwas_n_data[, c("PHENO", "N")]
    gwas_n_data <- setNames(gwas_n_data$N, gwas_n_data$PHENO)
  } else {
    if (!is.null(gwas_n_data)) {
      rlang::abort("Error: In binary or quantitative GWAS with known SDY, 'gwas_n_data' does not need to be provided")
    }
  }

  ##check eqtl_n_data
  if (!is.null(eqtl_n_data)) {
    if (!is.character(eqtl_n_data) || length(eqtl_n_data) != 1 || !file.exists(eqtl_n_data)) {
      rlang::abort("Error: 'eqtl_n_data' must be a non-null character string pointing to a readable file")
    }
    eqtl_n_data <- data.table::fread(eqtl_n_data, data.table=F, stringsAsFactors = F)
    if (!is.data.frame(eqtl_n_data) && !is.matrix(eqtl_n_data)) {
      rlang::abort("'eqtl_n_data' must be a valid file path or a data frame/matrix.")
    }
    colnames_upper <- toupper(colnames(eqtl_n_data))
    if (!all(c("TISSUE", "N") %in% colnames_upper)) {
      rlang::abort("Error: eqtl_n_data must contain columns 'TISSUE' and 'N' (case-insensitive).")
    }
    colnames(eqtl_n_data) <- colnames_upper
    eqtl_n_data <- eqtl_n_data[, c("TISSUE", "N")]
    eqtl_n_data <- setNames(eqtl_n_data$N, eqtl_n_data$TISSUE)
  } else {
    eqtl_n_data <- data.table::fread(system.file("extdata", "eqtl_samplesize.csv.gz", package = "mds2g"), data.table=F, stringsAsFactors = F)
    if (!is.data.frame(eqtl_n_data) && !is.matrix(eqtl_n_data)) {
      rlang::abort("'eqtl_n_data' must be a valid file path or a data frame/matrix.")
    }
    colnames_upper <- toupper(colnames(eqtl_n_data))
    if (!all(c("TISSUE", "N") %in% colnames_upper)) {
      rlang::abort("Error: eqtl_n_data must contain columns 'TISSUE' and 'N' (case-insensitive).")
    }
    colnames(eqtl_n_data) <- colnames_upper
    eqtl_n_data <- eqtl_n_data[, c("TISSUE", "N")]
    eqtl_n_data <- setNames(eqtl_n_data$N, eqtl_n_data$TISSUE)
  }

  ##check save_folder
  if (missing(save_folder) || !is.character(save_folder) || length(save_folder) != 1) {
    rlang::abort("'save_folder' must be provided and must be a valid directory path.")
  }

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(save_folder)) {
    dir.create(save_folder, recursive = TRUE)  # The `recursive = TRUE` allows the creation of nested directories.
  }

  # Check write permissions after ensuring the directory exists
  if (!file.access(save_folder, 2) == 0) {
    rlang::abort("'save_folder' must be a directory where the user has write permissions.")
  }

  save_folder <- normalizePath(save_folder)

  ##check work_part
  if (!is.null(work_part)) {
    if (!is.character(work_part) || length(work_part) != 1) {
      rlang::abort("Error: 'work_part' must be a single string in the format 'A/B'.")
    }
    if (!grepl("^\\d+/\\d+$", work_part)) {
      rlang::abort("Error: 'work_part' must be in the format 'A/B', where A and B are positive integers.")
    }
    work_parts <- as.numeric(unlist(strsplit(work_part, "/")))
    current_part <- work_parts[1]
    all_part <- work_parts[2]
    if (current_part <= 0 || all_part <= 0 || current_part > B) {
      rlang::abort("Error: 'A' and 'B' must be positive integers, and 'A' must be less than or equal to 'B'.")
    }
    splits <- split(SNP_data, sort(rep(1:all_part, length.out = nrow(SNP_data))))
    SNP_data <- splits[[current_part]]
  }

  ##check coloc_p
  if (!is.null(coloc_p1)) {
    if (!is_numeric(coloc_p1) && is.na(as.numeric(coloc_p1))) {
      stop("Error: 'coloc_p1' must be a numeric value.")
    }
    coloc_p1 <- as.numeric(coloc_p1)
  } else {
    coloc_p1 <- 1e-4
  }

  if (!is.null(coloc_p2)) {
    if (!is_numeric(coloc_p2) && is.na(as.numeric(coloc_p2))) {
      stop("Error: 'coloc_p2' must be a numeric value.")
    }
    coloc_p2 <- as.numeric(coloc_p2)
  } else {
    coloc_p2 <- 1e-4
  }

  if (!is.null(coloc_p12)) {
    if (!is_numeric(coloc_p12) && is.na(as.numeric(coloc_p12))) {
      stop("Error: 'coloc_p12' must be a numeric value.")
    }
    coloc_p12 <- as.numeric(coloc_p12)
  } else {
    coloc_p12 <- 1e-5
  }



  complement <- function(base) {
    return(chartr("ATGC", "TACG", base))
  }


  is_palindromic <- function(A1, A2) {
    return((A1 == complement(A2)) & (A2 == complement(A1)))
  }




  harmonize <- function(df) {
    rows_to_remove <- vector("logical", nrow(df))
    for (i in 1:nrow(df)) {
      if (is_palindromic(df$A1_gwas[i], df$A2_gwas[i]) | is_palindromic(df$A1_eqtl[i], df$A2_eqtl[i])) {
        rows_to_remove[i] <- TRUE
        next
      }

      if ((df$A1_gwas[i] == complement(df$A1_eqtl[i]) & df$A2_gwas[i] == complement(df$A2_eqtl[i])) |
          (df$A1_gwas[i] == complement(df$A2_eqtl[i]) & df$A2_gwas[i] == complement(df$A1_eqtl[i]))) {
        df$A1_eqtl[i] <- complement(df$A1_eqtl[i])
        df$A2_eqtl[i] <- complement(df$A2_eqtl[i])
      }

      if (df$A1_gwas[i] == df$A1_eqtl[i] & df$A2_gwas[i] == df$A2_eqtl[i]) {
        next
      } else if (df$A1_gwas[i] == df$A2_eqtl[i] & df$A2_gwas[i] == df$A1_eqtl[i]) {
        df$BETA_eqtl[i] <- -df$BETA_eqtl[i]
        temp <- df$A1_eqtl[i]
        df$A1_eqtl[i] <- df$A2_eqtl[i]
        df$A2_eqtl[i] <- temp
      } else {
        rows_to_remove[i] <- TRUE
      }
    }

    df <- df[!rows_to_remove, ]

    return(df)
  }





  ###coloc for six parameter combination

  if (gwas_data_type == "quant" && gwas_sdy_input == TRUE && gwas_beta_input == TRUE) {
    ##coloc
    coloc_merge_out <- data.frame()
    SNP_merge_out <- data.frame()
    required_cols <- c("SNP", "CHR", "BP", "BETA", "SE")
    if (strand_harmonize == TRUE) {
      required_cols <- c(required_cols, "A1", "A2")
    }

    for (tissue in tissue_list) {
      tissue_sample_size <- eqtl_n_data[[tissue]]

      for (i in 1:nrow(SNP_data)) {
        PHENO <- SNP_data$PHENO[i]
        SNP <- SNP_data$SNP[i]
        GENE_columns <- SNP_data[i, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6")]
        valid_genes <- GENE_columns[GENE_columns != "NONE"]

        gwas_sdy_data <- gwas_sdy_data[[PHENO]]

        gwas_file <- file.path(gwas_separate_data_folder, PHENO, paste0(SNP, ".txt"))
        if (!file.exists(gwas_file)) {
          warning(paste0(gwas_file, " does not exist, please check!"))
          next
        }
        gwas_data <- data.table::fread(gwas_file, sep="\t", data.table = F, stringsAsFactors = F)
        if (!is.data.frame(gwas_data) && !is.matrix(gwas_data)) {
          warning(paste0(gwas_file, " is not a valid file path point to a data frame/matrix, please check!"))
          next
        }
        if (!setequal(required_cols, colnames(gwas_data))) {
          missing_cols <- setdiff(required_cols, colnames(gwas_data))
          extra_cols <- setdiff(colnames(gwas_data), required_cols)
          error_message <- "Error: GWAS data columns do not match the required columns. "
          if (length(missing_cols) > 0) {
            error_message <- paste0(error_message, "The following required columns are missing: ", paste(missing_cols, collapse = ", "), ". ")
          }
          if (length(extra_cols) > 0) {
            error_message <- paste0(error_message, "The following columns are extra: ", paste(extra_cols, collapse = ", "), ". ")
          }
          rlang::abort(paste0(error_message, "Please confirm if the parameter settings about GWAS are the same as in the previous step!"))
        }

        dup_snp <- duplicated(gwas_data$SNP)
        gwas_data <- gwas_data[!dup_snp,]

        for (gene in valid_genes) {
          eqtl_file <- file.path(eqtl_separate_data_folder, tissue, "separate_file", paste0(tissue, "_", gene, ".txt"))

          eqtl_data <- tryCatch({
            data.table::fread(eqtl_file, sep="\t", header=TRUE, data.table = F, stringsAsFactors = F)
          }, error = function(e) {
            message(paste("Error in reading", eqtl_file, "- skipping this gene."))
            return(NULL)  # Return NULL if there is an error
          })

          if (is.null(eqtl_data)) next  # If eqtl_data is NULL, skip to the next iteration

          merged_df <- merge(gwas_data, eqtl_data, by="SNP", suffixes=c("_gwas", "_eqtl"))
          if (dim(merged_df)[1]==0) {
            message(paste("Skipping current iteration: The merged dataframe is empty. Please check if the SNP names are consistent between", gwas_file, "and", eqtl_file, "!"))
            next
          }
          dup_snp <- duplicated(merged_df$SNP)
          merged_df <- merged_df[!dup_snp,]
          merged_df <- merged_df[!is.na(merged_df$MAF), ]
          if (is.na(merged_df$MAF[1])) {
            message(paste0("Skipping current iteration: 'MAF' column is NA after match in ", eqtl_file))
            next
          }
          merged_df$varbeta_gwas <- merged_df$SE_gwas^2
          merged_df$varbeta_eqtl <- merged_df$SE_eqtl^2
          merged_df$MAF[merged_df$MAF == 0] <- 0.001
          if (strand_harmonize) {
            merged_df<-harmonize(merged_df)
          }

          coloc_gwas <- list(
            beta = merged_df$BETA_gwas,
            varbeta = merged_df$varbeta_gwas,
            type = "quant",
            sdY = gwas_sdy_data,
            snp = merged_df$SNP
          )

          coloc_eqtl <- list(
            beta = merged_df$BETA_eqtl,
            varbeta = merged_df$varbeta_eqtl,
            N = tissue_sample_size,
            type = "quant",
            MAF = merged_df$MAF,
            snp = merged_df$SNP
          )

          coloc.res <- coloc::coloc.abf(coloc_gwas, coloc_eqtl, p1=coloc_p1, p2=coloc_p2, p12=coloc_p12)

          coloc_summary <- coloc.res$summary
          coloc_summary$PHENO <- PHENO
          coloc_summary$SNP <- SNP
          coloc_summary$GENE <- gene
          coloc_summary$TISSUE <- tissue
          coloc_merge_out <- data.table::rbindlist(list(coloc_merge_out, coloc_summary), fill=TRUE)


          if (coloc.res$summary[6] >= 0.8) {
            SNP_summary <- subset(coloc.res$results, SNP.PP.H4 >= 0.8)
            if (nrow(SNP_summary) > 0) {
              SNP_summary$PHENO <- PHENO
              SNP_summary$SNP <- SNP
              SNP_summary$GENE <- gene
              SNP_summary$TISSUE <- tissue
              SNP_merge_out <- data.table::rbindlist(list(SNP_merge_out, SNP_summary), fill=TRUE)
            } else {
              warning("No rows found with SNP.PP.H4 >= 0.8")
              next
            }
          }
        }
      }
    }



  } else if (gwas_data_type == "quant" && gwas_sdy_input == FALSE && gwas_beta_input == TRUE) {
    coloc_merge_out <- data.frame()
    SNP_merge_out <- data.frame()
    required_cols <- c("SNP", "CHR", "BP", "BETA", "SE", "MAF")
    if (strand_harmonize == TRUE) {
      required_cols <- c(required_cols, "A1", "A2")
    }

    for (tissue in tissue_list) {
      tissue_sample_size <- eqtl_n_data[[tissue]]

      for (i in 1:nrow(SNP_data)) {
        PHENO <- SNP_data$PHENO[i]
        SNP <- SNP_data$SNP[i]
        GENE_columns <- SNP_data[i, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6")]
        valid_genes <- GENE_columns[GENE_columns != "NONE"]

        gwas_sample_size <- gwas_n_data[[PHENO]]

        gwas_file <- file.path(gwas_separate_data_folder, PHENO, paste0(SNP, ".txt"))
        if (!file.exists(gwas_file)) {
          warning(paste0(gwas_file, " does not exist, please check!"))
          next
        }
        gwas_data <- data.table::fread(gwas_file, sep="\t", data.table = F, stringsAsFactors = F)
        if (!is.data.frame(gwas_data) && !is.matrix(gwas_data)) {
          warning(paste0(gwas_file, " is not a valid file path point to a data frame/matrix, please check!"))
          next
        }
        if (!setequal(required_cols, colnames(gwas_data))) {
          missing_cols <- setdiff(required_cols, colnames(gwas_data))
          extra_cols <- setdiff(colnames(gwas_data), required_cols)
          error_message <- "Error: GWAS data columns do not match the required columns. "
          if (length(missing_cols) > 0) {
            error_message <- paste0(error_message, "The following required columns are missing: ", paste(missing_cols, collapse = ", "), ". ")
          }
          if (length(extra_cols) > 0) {
            error_message <- paste0(error_message, "The following columns are extra: ", paste(extra_cols, collapse = ", "), ". ")
          }
          rlang::abort(paste0(error_message, "Please confirm if the parameter settings about GWAS are the same as in the previous step!"))
        }

        dup_snp <- duplicated(gwas_data$SNP)
        gwas_data <- gwas_data[!dup_snp,]

        for (gene in valid_genes) {
          eqtl_file <- file.path(eqtl_separate_data_folder, tissue, "separate_file", paste0(tissue, "_", gene, ".txt"))

          eqtl_data <- tryCatch({
            data.table::fread(eqtl_file, sep="\t", header=TRUE, data.table = F, stringsAsFactors = F)
          }, error = function(e) {
            message(paste("Error in reading", eqtl_file, "- skipping this gene."))
            return(NULL)  # Return NULL if there is an error
          })

          if (is.null(eqtl_data)) next  # If eqtl_data is NULL, skip to the next iteration

          merged_df <- merge(gwas_data, eqtl_data, by="SNP", suffixes=c("_gwas", "_eqtl"))
          if (dim(merged_df)[1]==0) {
            message(paste("Skipping current iteration: The merged dataframe is empty. Please check if the SNP names are consistent between", gwas_file, "and", eqtl_file, "!"))
            next
          }
          dup_snp <- duplicated(merged_df$SNP)
          merged_df <- merged_df[!dup_snp,]
          merged_df <- merged_df[!is.na(merged_df$MAF_eqtl), ]
          if (is.na(merged_df$MAF_eqtl[1])) {
            message(paste0("Skipping current iteration: 'MAF_eqtl' column is NA after match in ", eqtl_file))
            next
          }
          merged_df$varbeta_gwas <- merged_df$SE_gwas^2
          merged_df$varbeta_eqtl <- merged_df$SE_eqtl^2
          merged_df$MAF_gwas[merged_df$MAF_gwas == 0] <- 0.001
          merged_df$MAF_eqtl[merged_df$MAF_eqtl == 0] <- 0.001
          if (strand_harmonize) {
            merged_df<-harmonize(merged_df)
          }

          coloc_gwas <- list(
            beta = merged_df$BETA_gwas,
            varbeta = merged_df$varbeta_gwas,
            N = gwas_sample_size,
            type = "quant",
            MAF = merged_df$MAF_gwas,
            snp = merged_df$SNP
          )

          coloc_eqtl <- list(
            beta = merged_df$BETA_eqtl,
            varbeta = merged_df$varbeta_eqtl,
            N = tissue_sample_size,
            type = "quant",
            MAF = merged_df$MAF_eqtl,
            snp = merged_df$SNP
          )

          coloc.res <- coloc::coloc.abf(coloc_gwas, coloc_eqtl, p1=coloc_p1, p2=coloc_p2, p12=coloc_p12)

          coloc_summary <- coloc.res$summary
          coloc_summary$PHENO <- PHENO
          coloc_summary$SNP <- SNP
          coloc_summary$GENE <- gene
          coloc_summary$TISSUE <- tissue
          coloc_merge_out <- data.table::rbindlist(list(coloc_merge_out, coloc_summary), fill=TRUE)


          if (coloc.res$summary[6] >= 0.8) {
            SNP_summary <- subset(coloc.res$results, SNP.PP.H4 >= 0.8)
            if (nrow(SNP_summary) > 0) {
              SNP_summary$PHENO <- PHENO
              SNP_summary$SNP <- SNP
              SNP_summary$GENE <- gene
              SNP_summary$TISSUE <- tissue
              SNP_merge_out <- data.table::rbindlist(list(SNP_merge_out, SNP_summary), fill=TRUE)
            } else {
              warning("No rows found with SNP.PP.H4 >= 0.8")
              next
            }
          }
        }
      }
    }





  } else if (gwas_data_type == "quant" && gwas_sdy_input == FALSE && gwas_beta_input == FALSE) {
    coloc_merge_out <- data.frame()
    SNP_merge_out <- data.frame()
    required_cols <- c("SNP", "CHR", "BP", "P", "MAF")
    if (strand_harmonize == TRUE) {
      required_cols <- c(required_cols, "A1", "A2")
    }

    for (tissue in tissue_list) {
      tissue_sample_size <- eqtl_n_data[[tissue]]

      for (i in 1:nrow(SNP_data)) {
        PHENO <- SNP_data$PHENO[i]
        SNP <- SNP_data$SNP[i]
        GENE_columns <- SNP_data[i, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6")]
        valid_genes <- GENE_columns[GENE_columns != "NONE"]

        gwas_sample_size <- gwas_n_data[[PHENO]]

        gwas_file <- file.path(gwas_separate_data_folder, PHENO, paste0(SNP, ".txt"))
        if (!file.exists(gwas_file)) {
          warning(paste0(gwas_file, " does not exist, please check!"))
          next
        }
        gwas_data <- data.table::fread(gwas_file, sep="\t", data.table = F, stringsAsFactors = F)
        if (!is.data.frame(gwas_data) && !is.matrix(gwas_data)) {
          warning(paste0(gwas_file, " is not a valid file path point to a data frame/matrix, please check!"))
          next
        }
        if (!setequal(required_cols, colnames(gwas_data))) {
          missing_cols <- setdiff(required_cols, colnames(gwas_data))
          extra_cols <- setdiff(colnames(gwas_data), required_cols)
          error_message <- "Error: GWAS data columns do not match the required columns. "
          if (length(missing_cols) > 0) {
            error_message <- paste0(error_message, "The following required columns are missing: ", paste(missing_cols, collapse = ", "), ". ")
          }
          if (length(extra_cols) > 0) {
            error_message <- paste0(error_message, "The following columns are extra: ", paste(extra_cols, collapse = ", "), ". ")
          }
          rlang::abort(paste0(error_message, "Please confirm if the parameter settings about GWAS are the same as in the previous step!"))
        }

        dup_snp <- duplicated(gwas_data$SNP)
        gwas_data <- gwas_data[!dup_snp,]

        for (gene in valid_genes) {
          eqtl_file <- file.path(eqtl_separate_data_folder, tissue, "separate_file", paste0(tissue, "_", gene, ".txt"))

          eqtl_data <- tryCatch({
            data.table::fread(eqtl_file, sep="\t", header=TRUE, data.table = F, stringsAsFactors = F)
          }, error = function(e) {
            message(paste("Error in reading", eqtl_file, "- skipping this gene."))
            return(NULL)  # Return NULL if there is an error
          })

          if (is.null(eqtl_data)) next  # If eqtl_data is NULL, skip to the next iteration

          merged_df <- merge(gwas_data, eqtl_data, by="SNP", suffixes=c("_gwas", "_eqtl"))
          if (dim(merged_df)[1]==0) {
            message(paste("Skipping current iteration: The merged dataframe is empty. Please check if the SNP names are consistent between", gwas_file, "and", eqtl_file, "!"))
            next
          }
          dup_snp <- duplicated(merged_df$SNP)
          merged_df <- merged_df[!dup_snp,]
          merged_df <- merged_df[!is.na(merged_df$MAF_eqtl), ]
          colnames(merged_df)[colnames(merged_df) == "BETA"] <- "BETA_eqtl"
          if (is.na(merged_df$MAF_eqtl[1])) {
            message(paste0("Skipping current iteration: 'MAF_eqtl' column is NA after match in ", eqtl_file))
            next
          }
          merged_df$varbeta_eqtl <- merged_df$SE^2
          merged_df$MAF_gwas[merged_df$MAF_gwas == 0] <- 0.001
          merged_df$MAF_eqtl[merged_df$MAF_eqtl == 0] <- 0.001
          if (strand_harmonize) {
            merged_df<-harmonize(merged_df)
          }

          coloc_gwas <- list(
            N = gwas_sample_size,
            type = "quant",
            MAF = merged_df$MAF_gwas,
            snp = merged_df$SNP,
            pvalues = merged_df$P_gwas
          )

          coloc_eqtl <- list(
            beta = merged_df$BETA_eqtl,
            varbeta = merged_df$varbeta_eqtl,
            N = tissue_sample_size,
            type = "quant",
            MAF = merged_df$MAF_eqtl,
            snp = merged_df$SNP
          )

          coloc.res <- coloc::coloc.abf(coloc_gwas, coloc_eqtl, p1=coloc_p1, p2=coloc_p2, p12=coloc_p12)

          coloc_summary <- coloc.res$summary
          coloc_summary$PHENO <- PHENO
          coloc_summary$SNP <- SNP
          coloc_summary$GENE <- gene
          coloc_summary$TISSUE <- tissue
          coloc_merge_out <- data.table::rbindlist(list(coloc_merge_out, coloc_summary), fill=TRUE)


          if (coloc.res$summary[6] >= 0.8) {
            SNP_summary <- subset(coloc.res$results, SNP.PP.H4 >= 0.8)
            if (nrow(SNP_summary) > 0) {
              SNP_summary$PHENO <- PHENO
              SNP_summary$SNP <- SNP
              SNP_summary$GENE <- gene
              SNP_summary$TISSUE <- tissue
              SNP_merge_out <- data.table::rbindlist(list(SNP_merge_out, SNP_summary), fill=TRUE)
            } else {
              warning("No rows found with SNP.PP.H4 >= 0.8")
              next
            }
          }
        }
      }
    }









  } else if (gwas_data_type == "quant" && gwas_sdy_input == TRUE && gwas_beta_input == FALSE) {
    coloc_merge_out <- data.frame()
    SNP_merge_out <- data.frame()
    required_cols <- c("SNP", "CHR", "BP", "P", "MAF")
    if (strand_harmonize == TRUE) {
      required_cols <- c(required_cols, "A1", "A2")
    }

    for (tissue in tissue_list) {
      tissue_sample_size <- eqtl_n_data[[tissue]]

      for (i in 1:nrow(SNP_data)) {
        PHENO <- SNP_data$PHENO[i]
        SNP <- SNP_data$SNP[i]
        GENE_columns <- SNP_data[i, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6")]
        valid_genes <- GENE_columns[GENE_columns != "NONE"]

        gwas_sdy_data <- gwas_sdy_data[[PHENO]]

        gwas_file <- file.path(gwas_separate_data_folder, PHENO, paste0(SNP, ".txt"))
        if (!file.exists(gwas_file)) {
          warning(paste0(gwas_file, " does not exist, please check!"))
          next
        }
        gwas_data <- data.table::fread(gwas_file, sep="\t", data.table = F, stringsAsFactors = F)
        if (!is.data.frame(gwas_data) && !is.matrix(gwas_data)) {
          warning(paste0(gwas_file, " is not a valid file path point to a data frame/matrix, please check!"))
          next
        }
        if (!setequal(required_cols, colnames(gwas_data))) {
          missing_cols <- setdiff(required_cols, colnames(gwas_data))
          extra_cols <- setdiff(colnames(gwas_data), required_cols)
          error_message <- "Error: GWAS data columns do not match the required columns. "
          if (length(missing_cols) > 0) {
            error_message <- paste0(error_message, "The following required columns are missing: ", paste(missing_cols, collapse = ", "), ". ")
          }
          if (length(extra_cols) > 0) {
            error_message <- paste0(error_message, "The following columns are extra: ", paste(extra_cols, collapse = ", "), ". ")
          }
          rlang::abort(paste0(error_message, "Please confirm if the parameter settings about GWAS are the same as in the previous step!"))
        }

        dup_snp <- duplicated(gwas_data$SNP)
        gwas_data <- gwas_data[!dup_snp,]

        for (gene in valid_genes) {
          eqtl_file <- file.path(eqtl_separate_data_folder, tissue, "separate_file", paste0(tissue, "_", gene, ".txt"))

          eqtl_data <- tryCatch({
            data.table::fread(eqtl_file, sep="\t", header=TRUE, data.table = F, stringsAsFactors = F)
          }, error = function(e) {
            message(paste("Error in reading", eqtl_file, "- skipping this gene."))
            return(NULL)  # Return NULL if there is an error
          })

          if (is.null(eqtl_data)) next  # If eqtl_data is NULL, skip to the next iteration

          merged_df <- merge(gwas_data, eqtl_data, by="SNP", suffixes=c("_gwas", "_eqtl"))
          if (dim(merged_df)[1]==0) {
            message(paste("Skipping current iteration: The merged dataframe is empty. Please check if the SNP names are consistent between", gwas_file, "and", eqtl_file, "!"))
            next
          }
          colnames(merged_df)[colnames(merged_df) == "BETA"] <- "BETA_eqtl"
          dup_snp <- duplicated(merged_df$SNP)
          merged_df <- merged_df[!dup_snp,]
          merged_df <- merged_df[!is.na(merged_df$MAF_eqtl), ]
          if (is.na(merged_df$MAF_eqtl[1])) {
            message(paste0("Skipping current iteration: 'MAF_eqtl' column is NA after match in ", eqtl_file))
            next
          }

          merged_df$varbeta_eqtl <- merged_df$SE^2
          merged_df$MAF_gwas[merged_df$MAF_gwas == 0] <- 0.001
          merged_df$MAF_eqtl[merged_df$MAF_eqtl == 0] <- 0.001
          if (strand_harmonize) {
            merged_df<-harmonize(merged_df)
          }

          coloc_gwas <- list(
            type = "quant",
            MAF = merged_df$MAF_gwas,
            snp = merged_df$SNP,
            sdy = gwas_sdy_data,
            pvalues = merged_df$P_gwas
          )

          coloc_eqtl <- list(
            beta = merged_df$BETA_eqtl,
            varbeta = merged_df$varbeta_eqtl,
            N = tissue_sample_size,
            type = "quant",
            MAF = merged_df$MAF_eqtl,
            snp = merged_df$SNP
          )

          coloc.res <- coloc::coloc.abf(coloc_gwas, coloc_eqtl, p1=coloc_p1, p2=coloc_p2, p12=coloc_p12)

          coloc_summary <- coloc.res$summary
          coloc_summary$PHENO <- PHENO
          coloc_summary$SNP <- SNP
          coloc_summary$GENE <- gene
          coloc_summary$TISSUE <- tissue
          coloc_merge_out <- data.table::rbindlist(list(coloc_merge_out, coloc_summary), fill=TRUE)


          if (coloc.res$summary[6] >= 0.8) {
            SNP_summary <- subset(coloc.res$results, SNP.PP.H4 >= 0.8)
            if (nrow(SNP_summary) > 0) {
              SNP_summary$PHENO <- PHENO
              SNP_summary$SNP <- SNP
              SNP_summary$GENE <- gene
              SNP_summary$TISSUE <- tissue
              SNP_merge_out <- data.table::rbindlist(list(SNP_merge_out, SNP_summary), fill=TRUE)
            } else {
              warning("No rows found with SNP.PP.H4 >= 0.8")
              next
            }
          }
        }
      }
    }

  } else if (gwas_data_type == "cc" && gwas_beta_input == TRUE) {
    coloc_merge_out <- data.frame()
    SNP_merge_out <- data.frame()
    required_cols <- c("SNP", "CHR", "BP", "BETA", "SE")
    if (strand_harmonize == TRUE) {
      required_cols <- c(required_cols, "A1", "A2")
    }

    for (tissue in tissue_list) {
      tissue_sample_size <- eqtl_n_data[[tissue]]

      for (i in 1:nrow(SNP_data)) {
        PHENO <- SNP_data$PHENO[i]
        SNP <- SNP_data$SNP[i]
        GENE_columns <- SNP_data[i, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6")]
        valid_genes <- GENE_columns[GENE_columns != "NONE"]

        gwas_s_data <- gwas_s_data[[PHENO]]

        gwas_file <- file.path(gwas_separate_data_folder, PHENO, paste0(SNP, ".txt"))
        if (!file.exists(gwas_file)) {
          warning(paste0(gwas_file, " does not exist, please check!"))
          next
        }
        gwas_data <- data.table::fread(gwas_file, sep="\t", data.table = F, stringsAsFactors = F)
        if (!is.data.frame(gwas_data) && !is.matrix(gwas_data)) {
          warning(paste0(gwas_file, " is not a valid file path point to a data frame/matrix, please check!"))
          next
        }
        if (!setequal(required_cols, colnames(gwas_data))) {
          missing_cols <- setdiff(required_cols, colnames(gwas_data))
          extra_cols <- setdiff(colnames(gwas_data), required_cols)
          error_message <- "Error: GWAS data columns do not match the required columns. "
          if (length(missing_cols) > 0) {
            error_message <- paste0(error_message, "The following required columns are missing: ", paste(missing_cols, collapse = ", "), ". ")
          }
          if (length(extra_cols) > 0) {
            error_message <- paste0(error_message, "The following columns are extra: ", paste(extra_cols, collapse = ", "), ". ")
          }
          rlang::abort(paste0(error_message, "Please confirm if the parameter settings about GWAS are the same as in the previous step!"))
        }

        dup_snp <- duplicated(gwas_data$SNP)
        gwas_data <- gwas_data[!dup_snp,]

        for (gene in valid_genes) {
          eqtl_file <- file.path(eqtl_separate_data_folder, tissue, "separate_file", paste0(tissue, "_", gene, ".txt"))

          eqtl_data <- tryCatch({
            data.table::fread(eqtl_file, sep="\t", header=TRUE, data.table = F, stringsAsFactors = F)
          }, error = function(e) {
            message(paste("Error in reading", eqtl_file, "- skipping this gene."))
            return(NULL)  # Return NULL if there is an error
          })

          if (is.null(eqtl_data)) next  # If eqtl_data is NULL, skip to the next iteration

          merged_df <- merge(gwas_data, eqtl_data, by="SNP", suffixes=c("_gwas", "_eqtl"))
          if (dim(merged_df)[1]==0) {
            message(paste("Skipping current iteration: The merged dataframe is empty. Please check if the SNP names are consistent between", gwas_file, "and", eqtl_file, "!"))
            next
          }
          dup_snp <- duplicated(merged_df$SNP)
          merged_df <- merged_df[!dup_snp,]
          merged_df <- merged_df[!is.na(merged_df$MAF), ]
          if (is.na(merged_df$MAF[1])) {
            message(paste0("Skipping current iteration: 'MAF_eqtl' column is NA after match in ", eqtl_file))
            next
          }
          merged_df$varbeta_gwas <- merged_df$SE_gwas^2
          merged_df$varbeta_eqtl <- merged_df$SE_eqtl^2
          merged_df$MAF[merged_df$MAF == 0] <- 0.001
          if (strand_harmonize) {
            merged_df<-harmonize(merged_df)
          }

          coloc_gwas <- list(
            beta = merged_df$BETA_gwas,
            varbeta = merged_df$varbeta_gwas,
            type = "cc",
            snp = merged_df$SNP,
            s = gwas_s_data
          )

          coloc_eqtl <- list(
            beta = merged_df$BETA_eqtl,
            varbeta = merged_df$varbeta_eqtl,
            N = tissue_sample_size,
            type = "quant",
            MAF = merged_df$MAF,
            snp = merged_df$SNP
          )

          coloc.res <- coloc::coloc.abf(coloc_gwas, coloc_eqtl, p1=coloc_p1, p2=coloc_p2, p12=coloc_p12)

          coloc_summary <- coloc.res$summary
          coloc_summary$PHENO <- PHENO
          coloc_summary$SNP <- SNP
          coloc_summary$GENE <- gene
          coloc_summary$TISSUE <- tissue
          coloc_merge_out <- data.table::rbindlist(list(coloc_merge_out, coloc_summary), fill=TRUE)


          if (coloc.res$summary[6] >= 0.8) {
            SNP_summary <- subset(coloc.res$results, SNP.PP.H4 >= 0.8)
            if (nrow(SNP_summary) > 0) {
              SNP_summary$PHENO <- PHENO
              SNP_summary$SNP <- SNP
              SNP_summary$GENE <- gene
              SNP_summary$TISSUE <- tissue
              SNP_merge_out <- data.table::rbindlist(list(SNP_merge_out, SNP_summary), fill=TRUE)
            } else {
              warning("No rows found with SNP.PP.H4 >= 0.8")
              next
            }
          }
        }
      }
    }

  } else if (gwas_data_type == "cc" && gwas_beta_input == FALSE) {
    coloc_merge_out <- data.frame()
    SNP_merge_out <- data.frame()
    required_cols <- c("SNP", "CHR", "BP", "P", "MAF")
    if (strand_harmonize == TRUE) {
      required_cols <- c(required_cols, "A1", "A2")
    }

    for (tissue in tissue_list) {
      tissue_sample_size <- eqtl_n_data[[tissue]]

      for (i in 1:nrow(SNP_data)) {
        PHENO <- SNP_data$PHENO[i]
        SNP <- SNP_data$SNP[i]
        GENE_columns <- SNP_data[i, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6")]
        valid_genes <- GENE_columns[GENE_columns != "NONE"]

        gwas_s_data <- gwas_s_data[[PHENO]]

        gwas_file <- file.path(gwas_separate_data_folder, PHENO, paste0(SNP, ".txt"))
        if (!file.exists(gwas_file)) {
          warning(paste0(gwas_file, " does not exist, please check!"))
          next
        }
        gwas_data <- data.table::fread(gwas_file, sep="\t", data.table = F, stringsAsFactors = F)
        if (!is.data.frame(gwas_data) && !is.matrix(gwas_data)) {
          warning(paste0(gwas_file, " is not a valid file path point to a data frame/matrix, please check!"))
          next
        }
        if (!setequal(required_cols, colnames(gwas_data))) {
          missing_cols <- setdiff(required_cols, colnames(gwas_data))
          extra_cols <- setdiff(colnames(gwas_data), required_cols)
          error_message <- "Error: GWAS data columns do not match the required columns. "
          if (length(missing_cols) > 0) {
            error_message <- paste0(error_message, "The following required columns are missing: ", paste(missing_cols, collapse = ", "), ". ")
          }
          if (length(extra_cols) > 0) {
            error_message <- paste0(error_message, "The following columns are extra: ", paste(extra_cols, collapse = ", "), ". ")
          }
          rlang::abort(paste0(error_message, "Please confirm if the parameter settings about GWAS are the same as in the previous step!"))
        }

        dup_snp <- duplicated(gwas_data$SNP)
        gwas_data <- gwas_data[!dup_snp,]

        for (gene in valid_genes) {
          eqtl_file <- file.path(eqtl_separate_data_folder, tissue, "separate_file", paste0(tissue, "_", gene, ".txt"))

          eqtl_data <- tryCatch({
            data.table::fread(eqtl_file, sep="\t", header=TRUE, data.table = F, stringsAsFactors = F)
          }, error = function(e) {
            message(paste("Error in reading", eqtl_file, "- skipping this gene."))
            return(NULL)  # Return NULL if there is an error
          })

          if (is.null(eqtl_data)) next  # If eqtl_data is NULL, skip to the next iteration

          merged_df <- merge(gwas_data, eqtl_data, by="SNP", suffixes=c("_gwas", "_eqtl"))
          if (dim(merged_df)[1]==0) {
            message(paste("Skipping current iteration: The merged dataframe is empty. Please check if the SNP names are consistent between", gwas_file, "and", eqtl_file, "!"))
            next
          }
          colnames(merged_df)[colnames(merged_df) == "BETA"] <- "BETA_eqtl"
          dup_snp <- duplicated(merged_df$SNP)
          merged_df <- merged_df[!dup_snp,]
          merged_df <- merged_df[!is.na(merged_df$MAF_eqtl), ]
          if (is.na(merged_df$MAF_eqtl[1])) {
            message(paste0("Skipping current iteration: 'MAF_eqtl' column is NA after match in ", eqtl_file))
            next
          }
          merged_df$varbeta_eqtl <- merged_df$SE^2
          merged_df$MAF_eqtl[merged_df$MAF_eqtl == 0] <- 0.001
          if (strand_harmonize) {
            merged_df<-harmonize(merged_df)
          }

          coloc_gwas <- list(
            MAF = merged_df$MAF_gwas,
            pvalues = merged_df$P_gwas,
            type = "cc",
            snp = merged_df$SNP,
            s = gwas_s_data
          )

          coloc_eqtl <- list(
            beta = merged_df$BETA_eqtl,
            varbeta = merged_df$varbeta_eqtl,
            N = tissue_sample_size,
            type = "quant",
            MAF = merged_df$MAF,
            snp = merged_df$SNP
          )

          coloc.res <- coloc::coloc.abf(coloc_gwas, coloc_eqtl, p1=coloc_p1, p2=coloc_p2, p12=coloc_p12)

          coloc_summary <- coloc.res$summary
          coloc_summary$PHENO <- PHENO
          coloc_summary$SNP <- SNP
          coloc_summary$GENE <- gene
          coloc_summary$TISSUE <- tissue
          coloc_merge_out <- data.table::rbindlist(list(coloc_merge_out, coloc_summary), fill=TRUE)


          if (coloc.res$summary[6] >= 0.8) {
            SNP_summary <- subset(coloc.res$results, SNP.PP.H4 >= 0.8)
            if (nrow(SNP_summary) > 0) {
              SNP_summary$PHENO <- PHENO
              SNP_summary$SNP <- SNP
              SNP_summary$GENE <- gene
              SNP_summary$TISSUE <- tissue
              SNP_merge_out <- data.table::rbindlist(list(SNP_merge_out, SNP_summary), fill=TRUE)
            } else {
              warning("No rows found with SNP.PP.H4 >= 0.8")
              next
            }
          }
        }
      }
    }
  } else {
    rlang::abort("Invalid combination of gwas_data_type, gwas_sdy_input, or gwas_beta_input.")
  }







  ###write result
  coloc_merge_out <- unique(coloc_merge_out)
  coloc_merge_out <- subset(coloc_merge_out, coloc_merge_out[[6]] >= 0.8)
  SNP_merge_out <- unique(SNP_merge_out)
  if(!is.null(work_part)) {
    output_file_coloc <- file.path(save_folder, paste0("coloc_part_", current_part, ".csv"))
    output_file_SNP <- file.path(save_folder, paste0("coloc_SNP_part_", current_part, ".csv"))
  } else {
    output_file_coloc <- file.path(save_folder, "coloc.csv")
    output_file_SNP <- file.path(save_folder, "coloc_SNP.csv")
  }
  write.csv(coloc_merge_out, output_file_coloc, row.names=F)
  write.csv(SNP_merge_out, output_file_SNP, row.names=F)
}












