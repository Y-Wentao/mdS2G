#' @title Clean the PoPs result
#'
#' @description
#' Use this function to convert gene names or delete specific genes (such as non-protein coding genes).
#'
#' @param pheno_name Character vector. Phenotypes to process.
#' @param SNP_data String (file path) or Data Frame. Must contain CHR, SNP, BP.
#' @param pops_result_folder String. Directory containing raw PoPS results.
#' @param pops_result_suffix String. Default ".preds".
#' @param save_folder String. Output directory.
#' @param annot_file String (path) or Data Frame. Gene annotation (CHR, START, END, GENE).
#' @param filter_gene We suggest deleting non coding genes to clarify the results. If you do not want to exclude, please select False. It's default to TRUE.
#' @param replace_list String (path) or Data Frame. Old ID -> New ID dictionary.
#'
#' @importFrom dplyr %>%
#' @importFrom data.table fread
#' @importFrom rlang abort warn inform
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom utils write.csv head
#' @importFrom stats quantile setNames
#' @export
pops_result_processing <- function(pheno_name,
                                   SNP_data,
                                   pops_result_folder,
                                   pops_result_suffix = NULL,
                                   save_folder,
                                   annot_file = NULL,
                                   filter_gene = TRUE,
                                   replace_list = NULL) {

  # --- 1. Input Validation & Preparation ---

  # Validate phenotype names
  pheno_list <- validate_pheno_name(pheno_name)

  rlang::inform(
    message = c(
      "i" = paste0("Phenotypes to process (", length(pheno_list), "):"),
      " " = paste(pheno_list, collapse = ", ")
    )
  )

  # Check directories
  if (missing(pops_result_folder) || !dir.exists(pops_result_folder)) {
    rlang::abort("'pops_result_folder' must be a valid directory path.")
  }
  if (missing(save_folder) || !dir.exists(save_folder)) {
    rlang::abort("'save_folder' must be a valid directory path.")
  }

  # Check suffix
  if (is.null(pops_result_suffix) || pops_result_suffix == "preds") {
    pops_result_suffix <- ".preds"
  }

  # Check filter_gene
  if (!is.logical(filter_gene) || length(filter_gene) != 1) {
    rlang::abort("'filter_gene' must be TRUE or FALSE.")
  }

  # --- 2. Process SNP Data ---

  # Note: causal_SNP = FALSE because PoPS logic typically relies on Position
  clean_snp_data <- process_snp_data(
    SNP_data = SNP_data,
    causal_SNP = FALSE,
    pheno_name = pheno_list
  )

  rlang::inform(paste("Total SNPs loaded and validated:", nrow(clean_snp_data)))

  # --- 3. Load Annotation File ---
  if (is.null(annot_file)) {
    # Default: Load internal package data
    annot_path <- system.file("extdata", "pops_gene_annot.csv.gz", package = "mds2g")
    if (annot_path == "") {
      rlang::abort("Default 'pops_gene_annot.csv.gz' not found in package 'mds2g'.")
    }
    annot_df <- data.table::fread(annot_path, data.table = FALSE)

  } else if (is.character(annot_file)) {
    # User file: existence check
    if (!file.exists(annot_file)) {
      rlang::abort(paste0("The provided annot_file does not exist: ", annot_file))
    }
    annot_df <- data.table::fread(annot_file, data.table = FALSE)

  } else {
    rlang::abort("'annot_file' must be NULL (default) or a valid file path string.")
  }

  # Standardize Columns
  req_annot_cols <- c("CHR", "START", "END", "SYMBOL")

  missing_cols <- setdiff(req_annot_cols, colnames(annot_df))
  if (length(missing_cols) > 0) {
    rlang::abort(c(
      "Invalid 'annot_file' format.",
      "x" = paste0("Missing required columns: ", paste(missing_cols, collapse = ", ")),
      "i" = "The file must contain: CHR, START, END, and SYMBOL."
    ))
  }

  annot_df <- annot_df[, req_annot_cols]
  rlang::inform(paste("Annotation loaded. Total genes:", nrow(annot_df)))


  # --- 4. load replace dictionary ---
  final_replace_list <- NULL

  if (!is.null(replace_list)) {
    if (is.character(replace_list)) {
      replace_df <- data.table::fread(replace_list, header = FALSE, data.table = FALSE)
    } else {
      replace_df <- replace_list
    }
  } else {
    r_path <- system.file("extdata", "pops_ID_replace_list.csv.gz", package = "mds2g")
    if (r_path != "") {
      replace_df <- data.table::fread(r_path, header = FALSE, data.table = FALSE)
    } else {
      replace_df <- NULL
    }
  }

  if (!is.null(replace_df)) {
    if (ncol(replace_df) < 2) rlang::abort("replace_list must have at least 2 columns.")
    final_replace_list <- setNames(replace_df[[2]], replace_df[[1]])
    rlang::inform(paste("Loaded ID replacement dictionary with", length(final_replace_list), "entries."))
  }


  # --- 5. Execution Loop ---

  for (pheno in pheno_list) {
    message(paste0("\n>>> Cleaning PoPS Result for: ", pheno))

    # Data Slicing
    current_SNP_data <- clean_snp_data[clean_snp_data$PHENO == pheno, ]

    if (nrow(current_SNP_data) == 0) {
      warning(paste("No SNP rows found for phenotype:", pheno, "- Skipping."))
      next
    }

    tryCatch({
      out_path <- process_pops(
        pheno = pheno,
        SNP_data = current_SNP_data,
        pops_folder = pops_result_folder,
        pops_result_suffix = pops_result_suffix,
        save_folder = save_folder,
        annot_file = annot_df,
        replace_list = final_replace_list,
        filter_gene = filter_gene
      )

      message(paste(">>> Success! Cleaned file saved to:", basename(out_path)))

    }, error = function(e) {
      message(paste("!!! Error processing", pheno, ":"))
      message(e$message)
    })
  }

  message("\nAll cleaning tasks completed.")
}
