#' Helper: Find Overlapping Genes within Window
#' Matches SNPs to Genes based on genomic coordinates.
#' @noRd
find_overlapping_genes <- function(snp_df, annot_df, window_size = 500000) {
  # Pre-calculate windows
  snp_df$START_WIN <- pmax(0, snp_df$BP - window_size)
  snp_df$END_WIN <- snp_df$BP + window_size
  snp_df$GENE_LIST <- NA_character_

  # Split by chromosome for performance
  snp_split <- split(snp_df, snp_df$CHR)
  annot_split <- split(annot_df, annot_df$CHR)

  results <- lapply(names(snp_split), function(chr) {
    curr_snps <- snp_split[[chr]]
    curr_annot <- annot_split[[chr]]

    if (is.null(curr_annot)) return(curr_snps)

    for (i in 1:nrow(curr_snps)) {
      s_start <- curr_snps$START_WIN[i]
      s_end <- curr_snps$END_WIN[i]

      # Find genes overlapping with the SNP window
      # annot_df now uses 'SYMBOL' column
      hits <- curr_annot[curr_annot$START <= s_end & curr_annot$END >= s_start, "SYMBOL"]

      if (length(hits) > 0) {
        curr_snps$GENE_LIST[i] <- paste(hits, collapse = ",")
      }
    }
    return(curr_snps)
  })

  if (length(results) == 0) return(snp_df)
  return(do.call(rbind, results))
}
