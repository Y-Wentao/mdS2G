#' Clean the PoPs result
#'
#' Use this function to convert gene names or delete  specific genes (such as non protein coding genes); You can specify the specific list to be converted or deleted by yourself
#'
#' @param pheno_name
#' @param SNP_data
#' @param pops_result_folder
#' @param save_folder
#' @param pops_result_suffix
#' @param annot_file
#' @param ID_replace
#' @param replace_list
#' @param ID_remove
#' @param remove_list
#'
#' @return
#' @export
#'
#' @examples
pops_result_processing <- function(pheno_name, SNP_data, pops_result_folder, save_folder, pops_result_suffix=NULL, annot_file = NULL,
                                   ID_replace = TRUE, replace_list = NULL, ID_remove = TRUE, remove_list =NULL)
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

  # check SNP_data
  check_and_standardize_SNP_data <- function(SNP_data, causal_SNP, pheno_name) {
    if (missing(SNP_data)) {
      rlang::abort("'SNP_data' must be provided as a file path.")
    }
    SNP_data <- data.table::fread(SNP_data, stringsAsFactors = FALSE, data.table=FALSE)
    if (!is.data.frame(SNP_data) && !is.matrix(SNP_data)) {
      rlang::abort("'SNP_data' must be a valid file path or a data frame/matrix.")
    }
    colnames_lower <- tolower(colnames(SNP_data))
    required_columns <- list(
      CHR = c("chr", "chrom"),
      SNP = "snp",
      BP = c("bp", "pos")
    )
    chr_col <- colnames(SNP_data)[colnames_lower %in% tolower(required_columns$CHR)]
    snp_col <- colnames(SNP_data)[colnames_lower %in% tolower(required_columns$SNP)]
    bp_col <- colnames(SNP_data)[colnames_lower %in% tolower(required_columns$BP)]
    if (length(chr_col) == 0 || length(snp_col) == 0 || length(bp_col) == 0) {
      rlang::abort("SNP_data must contain columns for 'CHR', 'SNP', and 'BP' (or equivalents).")
    }

    # check PHENO column when length of pheno_name > 1
    pheno_col <- NULL
    if (length(pheno_name) > 1) {
      pheno_candidates <- c("pheno", "phenotype", "trait")
      pheno_col <- colnames(SNP_data)[colnames_lower %in% pheno_candidates]
      if (length(pheno_col) == 0) {
        rlang::abort("SNP_data must contain a 'PHENO' (or equivalent) column when multiple 'pheno_name' are provided.")
      }

      selected_cols <- c(chr_col, snp_col, bp_col, pheno_col)
      SNP_data <- SNP_data[, selected_cols]
      colnames(SNP_data) <- c("CHR", "SNP", "BP", "PHENO")

      unique_pheno_in_data <- unique(SNP_data$PHENO)
      if (!all(unique(pheno_name) %in% unique_pheno_in_data)) {
        unmatched_pheno <- setdiff(pheno_name, unique_pheno_in_data)
        rlang::abort(paste("The following 'pheno_name' entries are not found in 'SNP_data':", paste(unmatched_pheno, collapse = ", ")))
      }
    } else {
      selected_cols <- c(chr_col, snp_col, bp_col)
      SNP_data <- SNP_data[, selected_cols]
      colnames(SNP_data) <- c("CHR", "SNP", "BP")
    }

    SNP_data$CHR <- gsub("chr", "", SNP_data$CHR, ignore.case = TRUE)  # remove 'chr' prefix
    SNP_data$CHR <- as.numeric(trimws(SNP_data$CHR))  # convert 'CHR' column to numeric
    SNP_data <- SNP_data[SNP_data$CHR %in% 1:22, ]  # only keep autosomal chromosome

    return(SNP_data)
  }

  SNP_data <- check_and_standardize_SNP_data(SNP_data, pheno_name = pheno_name)
  # check pops_result_folder
  if (missing(pops_result_folder) || !is.character(pops_result_folder) || length(pops_result_folder) != 1 || !dir.exists(pops_result_folder)) {
    rlang::abort("'pops_result_folder' must be provided and must be a valid directory path.")
  }
  pops_result_folder <- normalizePath(pops_result_folder)

  # check suffix
  # If pops_result_suffix is NULL or 'preds', assign it to '.preds'
  if (is.null(pops_result_suffix) || pops_result_suffix == "preds") {
    pops_result_suffix <- ".preds"
  } else if (!is.character(pops_result_suffix) || length(pops_result_suffix) != 1) {
    rlang::abort("'pops_result_suffix' must be a single string or NULL.")
  }


  # check save_folder
  if (missing(save_folder) || !is.character(save_folder) || length(save_folder) != 1 || !dir.exists(save_folder)) {
    rlang::abort("'save_folder' must be provided and must be a valid directory path.")
  }
  if (!file.access(save_folder, 2) == 0) {
    rlang::abort("'save_folder' must be a directory where the user has write permissions.")
  }
  save_folder <- normalizePath(save_folder)


  # Check ID_replace:
  if (!is.logical(ID_replace) || length(ID_replace) != 1) {
    rlang::abort("'ID_replace' must be either TRUE or FALSE.")
  }

  # Check ID_remove:
  if (!is.logical(ID_remove) || length(ID_remove) != 1) {
    rlang::abort("'ID_remove' must be either TRUE or FALSE.")
  }

  # check annot_file
  if (is.null(annot_file)) {
    annot_file <- data.table::fread(system.file("extdata", "GENE_annot.csv.gz", package = "mds2g"), header = TRUE, data.table = FALSE)
  } else {
    if(!is.character(annot_file) || length(annot_file) != 1) {
      rlang::abort("'annot_file' must be NULL (default) or a valid file.")
    }
    annot_file <- data.table::fread(annot_file, header = TRUE, data.table = FALSE)
    if (!(is.data.frame(annot_file) || is.matrix(annot_file))) {
      rlang::abort("'annot_file' must be a data frame or matrix.")
    }
  }

  colnames_lower <- tolower(colnames(annot_file))
  required_columns <- list(
    CHR = c("chr", "chrom"),
    START = "start",
    END = "end",
    GENE = c("ENSGID","GENE")
  )
  chr_col <- colnames(annot_file)[colnames_lower %in% tolower(required_columns$CHR)]
  start_col <- colnames(annot_file)[colnames_lower %in% tolower(required_columns$START)]
  end_col <- colnames(annot_file)[colnames_lower %in% tolower(required_columns$END)]
  gene_col <- colnames(annot_file)[colnames_lower %in% tolower(required_columns$GENE)]
  if (length(chr_col) == 0 || length(start_col) == 0 || length(end_col) == 0 || length(gene_col) == 0) {
    rlang::abort("annot_file must contain columns for 'CHR', 'START', 'END' and 'GENE'(or equivalents).")
  }
  selected_cols <- c(chr_col, start_col, end_col, gene_col)
  annot_file <- annot_file[, selected_cols]
  colnames(annot_file) <- c("CHR", "START", "END", "GENE")
  non_ensg_genes <- annot_file$GENE[!grepl("^ENSG", annot_file$GENE)]
  non_ensg_count <- length(non_ensg_genes)
  total_genes_count <- nrow(annot_file)
  if (non_ensg_count > 0) {
    warning(paste0(
      "In the provided annot_file, there are ", non_ensg_count,
      " out of ", total_genes_count,
      " gene names that do not conform to the ENSG ID coding rule. Please check!"
    ))
  }

  # check replace_list
  if (is.null(replace_list)) {
    replace_list <- data.table::fread(system.file("extdata", "pops_ID_replace_list.csv.gz", package = "mds2g"),header=F,data.table = F,stringsAsFactors = F)
  } else {
    if (!is.character(replace_list) || length(replace_list) != 1) {
      rlang::abort("'replace_list' must be NULL (default) or a valid file.")
    }
    replace_list <- data.table::fread(replace_list, header = FALSE, data.table = FALSE, stringsAsFactors = F)
    if (!(is.data.frame(replace_list) || is.matrix(replace_list))) {
      rlang::abort("'replace_list' must be a data frame or matrix.")
    }
    if (ncol(replace_list) != 2) {
      rlang::abort("'replace_list' must be a data frame or matrix with 2 columns (refer to old ID and new ID).")
    }
  }
  replace_list <- setNames(replace_list[[2]], replace_list[[1]])

  # check remove_list
  if (is.null(remove_list)) {
    remove_list <- data.table::fread(system.file("extdata", "pops_ID_remove_list.csv.gz", package = "mds2g"),header=F,data.table = F,stringsAsFactors = F)
  } else {
    if (!is.character(remove_list) || length(remove_list) != 1) {
      rlang::abort("'remove_list' must be NULL (default) or a valid file.")
    }
    remove_list <- data.table::fread(remove_list, header = FALSE, data.table = FALSE)

    # check remove_list
    if (!(is.data.frame(remove_list) || is.matrix(remove_list))) {
      rlang::abort("'remove_list' must be a data frame or matrix.")
    }
    if (ncol(replace_list) != 1) {
      rlang::abort("'replace_list' must be a data frame or matrix with single column")
    }
  }
  remove_list <- remove_list %>% dplyr::pull(1)

  result_processing <- function(pheno_name, pops_data, SNP_data, pops_result_folder, save_folder, pops_result_suffix = pops_result_suffix, annot_file = annot_file,
                                ID_replace = TRUE, replace_list = replace_list, ID_remove = TRUE, remove_list =remove_list)
  {
    pops_data <- pops_data %>%
      dplyr::mutate(ENSGID = recode(ENSGID, !!!replace_list))
    pops_data <- pops_data[!(pops_data$ENSGID %in% remove_list), ]
    score_thre <- sort(pops_data$PoPS_Score, decreasing = TRUE)[ceiling(0.20 * nrow(pops_data))]

    find_genes <- function(SNP_data, annot_file) {
      SNP_data$GENE <- NA
      SNP_data$START <- NA
      SNP_data$END <- NA
      for (i in 1:nrow(SNP_data)) {
        start_window <- max(0, SNP_data$BP[i] - 500000)
        end_window <- SNP_data$BP[i] + 500000
        SNP_data$START[i] <- start_window
        SNP_data$END[i] <- end_window
        overlapping_genes <- annot_file[
          annot_file$CHR == SNP_data$CHR[i] &
            ((annot_file$START >= start_window & annot_file$START <= end_window) |
               (annot_file$END >= start_window & annot_file$END <= end_window) |
               (annot_file$START <= start_window & annot_file$END >= end_window)),
          'GENE']
        SNP_data$GENE[i] <- paste(overlapping_genes, collapse = ',')
      }
      return(SNP_data)
    }
    SNP_data_with_genes <- find_genes(SNP_data, annot_file)
    SNP_data_with_genes$GENE1 <- NA
    SNP_data_with_genes$GENE2 <- NA
    SNP_data_with_genes$GENE3 <- NA
    SNP_data_with_genes$GENE4 <- NA
    SNP_data_with_genes$GENE5 <- NA
    for (i in 1:nrow(SNP_data_with_genes)) {
      genes <- unlist(strsplit(SNP_data_with_genes$GENE[i], ","))
      if (length(genes) == 0) {
        next
      }
      # get score
      gene_scores <- pops_data[pops_data$ENSGID %in% genes, ]
      gene_scores <- gene_scores[!is.na(gene_scores$PoPS_Score), ]
      # choose the top five gene
      top_genes <- gene_scores[order(-gene_scores$PoPS_Score), ][1:min(5, nrow(gene_scores)), ]
      # keep the gene which score higher than thres
      for (j in 1:5) {
        if (j <= nrow(top_genes) && top_genes$PoPS_Score[j] >= score_thre) {
          SNP_data_with_genes[i, paste0("GENE", j)] <- top_genes$ENSGID[j]
        } else {
          SNP_data_with_genes[i, paste0("GENE", j)] <- "NONE"
        }
      }
    }
    SNP_data_with_genes$PHENO <- NULL
    SNP_data_with_genes$`SNP_data$PHENO` <- NULL
    write.csv(SNP_data_with_genes,file.path(save_folder, paste0(pheno_name, "_pops_result_processed.csv")), row.names=F)
  }


  if (length(pheno_name) > 1) {
    for (pheno in pheno_name) {
      current_SNP_data <- subset(SNP_data, SNP_data$PHENO == pheno)
      current_pheno_name <- unique(current_SNP_data$PHENO)

      ##check pops_data
        file_path <- file.path(pops_result_folder, paste0(current_pheno_name, pops_result_suffix))
      if (!file.exists(file_path)) {
        rlang::abort(paste0("File does not exist: ", file_path))
      }
      pops_data <- data.table::fread(file_path, header = TRUE, stringsAsFactors = FALSE)
      if (!(is.data.frame(pops_data) || is.matrix(pops_data))) {
        rlang::abort(paste0("The file '", file_path, "' must be a data frame or matrix."))
      }
      expected_colnames <- c("ENSGID", "PoPS_Score", "Y", "Y_proj", "project_out_covariates_gene", "feature_selection_gene", "training_gene")
      actual_colnames <- colnames(pops_data)
      if (!identical(actual_colnames, expected_colnames)) {
        warning(paste0(
          "The expected pops_result file column names should be:\n",
          paste(expected_colnames, collapse = "\t"), "\n",
          "But the provided file column names are:\n",
          paste(actual_colnames, collapse = "\t"), "\n",
          "Please check the input file for any issues."
        ))
      }

      result_processing(
        pheno_name = current_pheno_name,
        pops_data = pops_data,
        SNP_data = current_SNP_data,
        pops_result_folder = pops_result_folder,
        save_folder = save_folder,
        pops_result_suffix = pops_result_suffix,
        annot_file = annot_file,
        ID_replace = TRUE,
        replace_list = replace_list,
        ID_remove = TRUE,
        remove_list = remove_list
      )
    }
  } else {
    file_path <- file.path(pops_result_folder, paste0(pheno_name, pops_result_suffix))
    if (!file.exists(file_path)) {
      rlang::abort(paste0("File does not exist: ", file_path))
    }
    pops_data <- data.table::fread(file_path, header = TRUE, stringsAsFactors = FALSE)
    if (!(is.data.frame(pops_data) || is.matrix(pops_data))) {
      rlang::abort(paste0("The file '", file_path, "' must be a data frame or matrix."))
    }
    expected_colnames <- c("ENSGID", "PoPS_Score", "Y", "Y_proj", "project_out_covariates_gene", "feature_selection_gene", "training_gene")
    actual_colnames <- colnames(pops_data)
    if (!identical(actual_colnames, expected_colnames)) {
      warning(paste0(
        "The expected pops_result file column names should be:\n",
        paste(expected_colnames, collapse = "\t"), "\n",
        "But the provided file column names are:\n",
        paste(actual_colnames, collapse = "\t"), "\n",
        "Please check the input file for any issues."
      ))
    }



    result_processing(
      pheno_name = pheno_name,
      pops_data = pops_data,
      SNP_data = SNP_data,
      pops_result_folder = pops_result_folder,
      save_folder = save_folder,
      pops_result_suffix = pops_result_suffix,
      annot_file = annot_file,
      ID_replace = TRUE,
      replace_list = replace_list,
      ID_remove = TRUE,
      remove_list = remove_list
    )
  }
}






