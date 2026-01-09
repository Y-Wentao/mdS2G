#' Extract nominated genes for colocalization-based linking
#'
#' Use this function to organize the list of potentially relevant genes obtained from locus-based linking (cS2G) and similarity-based linking (PoPs).
#' Low confidence cS2G genes (< 0.8) and all PoPs genes are nominated for colocalization.
#' High confidence cS2G genes (>= 0.8) are saved separately as credible genes.
#'
#' @param pheno_name Character string or vector. Comma-separated or vector of phenotype names.
#' @param cs2g_result_folder Path to the folder containing cS2G results.
#' @param pops_result_folder Path to the folder containing PoPs results.
#' @param save_folder Path where the output files will be saved.
#'
#' @export
nominated_gene_extract <- function(pheno_name, cs2g_result_folder, pops_result_folder, save_folder) {

  # 1. Validate and Standardize Inputs
  pheno_name <- validate_pheno_name(pheno_name)
  rlang::inform(
    message = c(
      "i" = paste0("Phenotypes to be analyzed (", length(pheno_name), "):"),
      " " = paste(pheno_name, collapse = ", ")
    )
  )

  # Check Directories
  cs2g_result_folder <- check_dir(cs2g_result_folder, "cs2g_result_folder")
  pops_result_folder <- check_dir(pops_result_folder, "pops_result_folder")
  save_folder        <- check_dir(save_folder, "save_folder")

  if (file.access(save_folder, 2) != 0) {
    rlang::abort("'save_folder' is not writable.")
  }

  # 2. Inner Processing Function
  process_single_pheno <- function(pheno) {

    # --- Load and Clean cS2G Data ---
    cs2g_file <- file.path(cs2g_result_folder, paste0(pheno, "_cS2G.csv"))
    if (!file.exists(cs2g_file)) rlang::abort(paste("cS2G file not found for:", pheno))

    cs2g_data <- data.table::fread(cs2g_file, data.table = FALSE, stringsAsFactors = FALSE)

    req_cols_cs2g <- c("SNP", "SYMBOL", "ENSGID", "RAW_SCORE", "WEIGHTED_SCORE", "STRATEGY", "PRECISION", "PIP", "CONFIDENCE_SCORE")
    missing_cs2g <- setdiff(req_cols_cs2g, colnames(cs2g_data))
    if (length(missing_cs2g) > 0) {
      rlang::abort(paste0("cS2G file for ", pheno, " missing columns: ", paste(missing_cs2g, collapse = ", ")))
    }
    cs2g_data <- cs2g_data[, req_cols_cs2g]

    # --- Load and Clean PoPs Data ---
    pops_file <- file.path(pops_result_folder, paste0(pheno, "_pops_result_processed.csv"))
    if (!file.exists(pops_file)) rlang::abort(paste("PoPs file not found for:", pheno))

    pops_data <- data.table::fread(pops_file, data.table = FALSE, stringsAsFactors = FALSE)

    req_cols_pops <- c("CHR", "SNP", "BP", "START", "END", paste0("GENE", 1:5)) # Assuming GENE1-5 exists
    missing_pops <- setdiff(req_cols_pops, colnames(pops_data))
    if (length(missing_pops) > 0) {
      rlang::abort(paste0("PoPs file for ", pheno, " missing columns: ", paste(missing_pops, collapse = ", ")))
    }
    pops_data <- pops_data[, req_cols_pops]

    # --- Logic: Separate High vs Low Confidence ---
    # Low confidence (< 0.8): To be merged with PoPs for colocalization
    cs2g_low_conf <- cs2g_data[cs2g_data$CONFIDENCE_SCORE < 0.8, ]

    # High confidence (>= 0.8): Credible genes, saved separately
    cs2g_high_conf <- cs2g_data[cs2g_data$CONFIDENCE_SCORE >= 0.8, ]
    cs2g_high_conf$PHENO <- pheno

    # --- Logic: Extract Gene List ---
    # Genes from Low Conf cS2G + All PoPs Genes
    pops_genes_vec <- unlist(pops_data[, paste0("GENE", 1:5)])

    candidate_genes <- unique(na.omit(c(
      cs2g_low_conf$ENSGID,
      pops_genes_vec
    )))
    # Remove "NONE" or empty strings
    candidate_genes <- candidate_genes[candidate_genes != "NONE" & candidate_genes != ""]

    # --- Logic: Merge Data for Output ---
    # Merge PoPs with Low Confidence cS2G info based on SNP
    merged_data <- merge(pops_data, cs2g_low_conf[, c("SNP", "ENSGID")], by = "SNP", all.x = TRUE)

    # Handle the merged cS2G ID (rename to GENE6)
    merged_data$GENE6 <- merged_data$ENSGID
    merged_data$GENE6[is.na(merged_data$GENE6)] <- "NONE"
    merged_data$ENSGID <- NULL # Remove the original column
    merged_data$PHENO <- pheno

    return(list(
      gene_list = candidate_genes,
      merged_data = merged_data,
      cs2g_cred = cs2g_high_conf
    ))
  }

  # 3. Execute Processing (Vectorized)
  message("Processing ", length(pheno_name), " phenotypes...")

  results_list <- lapply(pheno_name, process_single_pheno)

  # 4. Aggregate Results
  # Combine gene lists (and dedup)
  all_genes <- unique(unlist(lapply(results_list, function(x) x$gene_list)))

  # Combine DataFrames
  all_merged_data <- dplyr::bind_rows(lapply(results_list, function(x) x$merged_data))
  all_cs2g_cred   <- dplyr::bind_rows(lapply(results_list, function(x) x$cs2g_cred))

  # 6. Save Outputs
  # Save Gene List
  out_gene_file <- file.path(save_folder, "Gene_for_coloc.txt")
  write.table(all_genes, file = out_gene_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

  # Save Merged SNP Data
  out_snp_file <- file.path(save_folder, "SNP_data_for_coloc.csv")
  write.csv(all_merged_data, file = out_snp_file, row.names = FALSE)

  # Save Credible Genes
  out_cred_file <- file.path(save_folder, "cS2G_cred_gene.csv")
  write.csv(all_cs2g_cred, file = out_cred_file, row.names = FALSE)

  message("Processing completed.")
  message("  - Gene list saved to: ", out_gene_file)
  message("  - SNP data saved to:  ", out_snp_file)
  message("  - Credible genes saved to: ", out_cred_file)
}
