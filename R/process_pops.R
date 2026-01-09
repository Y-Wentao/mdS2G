#' Process PoPS Results for a Single Phenotype
#' Cleans PoPS IDs, filters genes, maps to SNPs, and selects top candidates.
#' @importFrom dplyr %>%
#' @noRd
process_pops <- function(pheno, SNP_data, pops_folder, pops_result_suffix, save_folder, annot_file, replace_list, filter_gene) {
  # 1. Read PoPS Result
  pops_file <- file.path(pops_folder, paste0(pheno, pops_result_suffix))
  if (!file.exists(pops_file)) rlang::abort(paste("PoPS result file missing:", pops_file))

  pops_data <- data.table::fread(pops_file, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)

  req_pops_cols <- c("ENSGID", "PoPS_Score")
  if (!all(req_pops_cols %in% colnames(pops_data))) {
    rlang::abort(paste("PoPS file must contain columns:", paste(req_pops_cols, collapse=", ")))
  }

  # 2. ID Cleaning & Dictionary Replacement
  if (!is.null(replace_list) && length(replace_list) > 0) {
    idx <- pops_data$ENSGID %in% names(replace_list)
    if (any(idx)) {
      pops_data$ENSGID[idx] <- replace_list[pops_data$ENSGID[idx]]
    }
  }

  unique_genes <- unique(pops_data$ENSGID)
  total_input_genes <- length(unique_genes)


  # 3. Query SYMBOL and GENETYPE using the standardized Ensembl IDs
  gene_info_df <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique_genes,
    columns = c("SYMBOL", "GENETYPE"),
    keytype = "ENSEMBL"
  ))
  gene_info_df <- gene_info_df[!duplicated(gene_info_df$ENSEMBL), ]

  # Remove unmatched, duplicated and non-protein-coding gene
  gene_info_df$REMOVAL_REASON <- NA_character_
  is_unmapped <- is.na(gene_info_df$SYMBOL)
  gene_info_df$REMOVAL_REASON[is_unmapped] <- "Unmapped_in_OrgDb"
  mask_mapped <- !is.na(gene_info_df$SYMBOL)
  is_dup_symbol <- duplicated(gene_info_df$SYMBOL) & mask_mapped
  to_mark_dup <- is_dup_symbol & is.na(gene_info_df$REMOVAL_REASON)
  gene_info_df$REMOVAL_REASON[to_mark_dup] <- "Duplicate_Symbol"

  if (filter_gene) {
    is_non_coding <- (is.na(gene_info_df$GENETYPE) | gene_info_df$GENETYPE != "protein-coding")
    to_mark_nc <- is_non_coding & is.na(gene_info_df$REMOVAL_REASON)
    gene_info_df$REMOVAL_REASON[to_mark_nc] <- "Non_Protein_Coding"
  }

  removed_df <- gene_info_df[!is.na(gene_info_df$REMOVAL_REASON), ]
  valid_map <- gene_info_df[is.na(gene_info_df$REMOVAL_REASON), ]

  # Summary
  n_unmapped <- sum(gene_info_df$REMOVAL_REASON == "Unmapped_in_OrgDb", na.rm = TRUE)
  n_dups <- sum(gene_info_df$REMOVAL_REASON == "Duplicate_Symbol", na.rm = TRUE)
  n_non_coding <- sum(gene_info_df$REMOVAL_REASON == "Non_Protein_Coding", na.rm = TRUE)

  n_total <- nrow(gene_info_df)
  n_mapped_unique <- n_total - n_unmapped - n_dups

  rlang::inform(c(
    "i" = paste0("Gene Filtering Summary for ", pheno, ":"),
    " " = paste0("  - Total unique genes: ", n_total),
    " " = paste0("  - Removed (Unmapped): ", n_unmapped, " (", round(n_unmapped/n_total*100, 2), "%)"),
    " " = paste0("  - Removed (Duplicate Symbol): ", n_dups),
    if(filter_gene) " " = paste0("  - Removed (Non-coding): ", n_non_coding, " (", round(n_non_coding/n_mapped_unique*100, 3), "% of mapped)"),
    "v" = paste0("  - Final Retained Genes: ", nrow(valid_map))
  ))

  # 4. Export Removed Genes
  if (nrow(removed_df) > 0) {
    out_removed <- removed_df[, c("ENSEMBL", "SYMBOL", "GENETYPE", "REMOVAL_REASON")]

    removed_file <- file.path(save_folder, paste0(pheno, "_pops_removed_genes.csv"))
    write.csv(out_removed, removed_file, row.names = FALSE)
    rlang::inform(paste0("  i Details of removed genes saved to: ", basename(removed_file)))
  }


  # 5. SNP Mapping & Top 5 Selection
  pops_cleaned <- merge(pops_data, valid_map, by.x = "ENSGID", by.y = "ENSEMBL")
  pops_cleaned <- pops_cleaned[, c("ENSGID", "SYMBOL", "PoPS_Score")]

  if (nrow(pops_cleaned) == 0) rlang::abort(paste("No genes left after filtering for pheno:", pheno))

  score_threshold <- quantile(pops_cleaned$PoPS_Score, 0.80, na.rm = TRUE)

  # Map SNPs to Genes using the provided annotation file (contains SYMBOL)
  snp_mapped <- find_overlapping_genes(SNP_data, annot_file)

  gene_cols <- paste0("GENE", 1:5)
  symbol_cols <- paste0("SYMBOL", 1:5)
  snp_mapped[c(gene_cols, symbol_cols)] <- "NONE"

  # Create lookup: SYMBOL -> Score
  gene_score_map <- setNames(pops_cleaned$PoPS_Score, pops_cleaned$SYMBOL)
  gene_id_map <- setNames(pops_cleaned$ENSGID, pops_cleaned$SYMBOL)

  for (i in 1:nrow(snp_mapped)) {
    genes_str <- snp_mapped$GENE_LIST[i]
    if (is.na(genes_str) || genes_str == "") next

    genes <- unlist(strsplit(genes_str, ","))

    # Retrieve scores for mapped symbols
    scores <- gene_score_map[genes]
    scores <- scores[!is.na(scores)]

    if (length(scores) == 0) next

    # Sort and select Top 5
    scores <- sort(scores, decreasing = TRUE)
    top_candidates <- head(scores, 5)

    for (j in seq_along(top_candidates)) {
      if (top_candidates[j] >= score_threshold) {
        current_symbol <- names(top_candidates)[j]
        current_ensgid <- gene_id_map[current_symbol]

        snp_mapped[i, gene_cols[j]] <- current_ensgid
        snp_mapped[i, symbol_cols[j]] <- current_symbol
      }
    }
  }

  # 6. Save Result
  final_cols <- c("CHR", "SNP", "BP", "START_WIN", "END_WIN", gene_cols, symbol_cols)
  snp_final <- snp_mapped[, final_cols, drop = FALSE] %>%
    dplyr::rename(START = START_WIN,
                  END = END_WIN)

  out_file <- file.path(save_folder, paste0(pheno, "_pops_result_processed.csv"))
  write.csv(snp_final, out_file, row.names = FALSE)

  return(out_file)
}
