#' Harmonize strands between GWAS and eQTL
#' @noRd
harmonize_data <- function(df) {
  # Columns required: A1_gwas, A2_gwas, A1_eqtl, A2_eqtl, BETA_eqtl

  # 1. Identify Palindromic (A/T, G/C) -> Remove
  is_pal <- function(a1, a2) {
    (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
      (a1 == "G" & a2 == "C") | (a1 == "C" & a2 == "G")
  }

  pal_mask <- is_pal(df$A1_gwas, df$A2_gwas) | is_pal(df$A1_eqtl, df$A2_eqtl)
  df <- df[!pal_mask, ]
  if(nrow(df) == 0) return(df)

  # 2. Complement Helper
  complement <- function(seq) {
    chartr("ATGC", "TACG", seq)
  }

  # 3. Check Alignment
  # Case A: Perfect match
  match_mask <- (df$A1_gwas == df$A1_eqtl) & (df$A2_gwas == df$A2_eqtl)

  # Case B: Swapped (Flip Beta)
  swap_mask <- (df$A1_gwas == df$A2_eqtl) & (df$A2_gwas == df$A1_eqtl)

  # Case C: Complement (Flip Strand)
  comp_mask <- (df$A1_gwas == complement(df$A1_eqtl)) & (df$A2_gwas == complement(df$A2_eqtl))

  # Case D: Swapped Complement (Flip Strand + Flip Beta)
  swap_comp_mask <- (df$A1_gwas == complement(df$A2_eqtl)) & (df$A2_gwas == complement(df$A1_eqtl))

  # Apply Changes
  # Flip Beta for swaps
  to_flip_beta <- swap_mask | swap_comp_mask
  df$BETA_eqtl[to_flip_beta] <- -df$BETA_eqtl[to_flip_beta]

  # For output correctness (optional, but good for debugging), update alleles
  # (Though coloc only cares about Beta matching the reference, updating alleles confirms alignment)
  df$A1_eqtl[to_flip_beta] <- df$A1_gwas[to_flip_beta]
  df$A2_eqtl[to_flip_beta] <- df$A2_gwas[to_flip_beta]

  # Keep only harmonizable rows
  keep_mask <- match_mask | swap_mask | comp_mask | swap_comp_mask
  df <- df[keep_mask, ]

  return(df)
}
