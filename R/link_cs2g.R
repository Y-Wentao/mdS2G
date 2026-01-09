#' Link SNPs to Genes using cS2G strategy (Internal)
#' @importFrom dplyr %>%
#' @noRd


link_cs2g <- function (pheno_name, SNP_data, s2g_folder, save_folder, keep_intermediate_file = FALSE,
                       intermediate_file_folder =FALSE, filter_gene = TRUE, custom_weights = FALSE,
                       custom_precisions = FALSE, causal_SNP = TRUE, assumed_PIP = NULL, keep_unmatched_gene = FALSE)
{
  pheno <- pheno_name
  dir_list <- c("Exon", "Promoter", "finemappedciseQTLs_eQTLGen",
                "finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood")

  #######################Step1: map SNP to GENE#######################################

  # Loop through each strategy
  raw_results_list <- lapply(dir_list, function(dir) {

    # --- 1. Read and Match SNPs ---

    # Split SNP data by chromosome to process in chunks
    # Assumes SNP_data columns: 1=CHR, 2=SNP, 3=POS
    snp_chunks <- split(SNP_data, SNP_data[[1]])

    strategy_results_list <- lapply(snp_chunks, function(chunk) {
      current_chr <- chunk[[1]][1]
      chr_gz_file <- normalizePath(file.path(s2g_folder, dir, paste0("chr", current_chr, ".bed.gz")))

      # Strict error checking for file existence
      if (!file.exists(chr_gz_file)) {
        rlang::abort(paste0("Annotation file missing: ", chr_gz_file))
      }

      # Read annotation file (Start, End, Gene, Value)
      chr_data <- data.table::fread(chr_gz_file, header = FALSE, sep = "\t",
                                    select = c(2, 3, 4, 5),
                                    col.names = c("start", "end", "gene", "value"),
                                    data.table = FALSE)

      # Find SNPs located within gene intervals
      result_rows <- list()
      for(k in 1:nrow(chunk)) {
        pos <- chunk[k, 3]    # SNP Position
        snp_id <- chunk[k, 2] # SNP ID
        hits <- chr_data[pos >= chr_data$start & pos <= chr_data$end, ]

        if (nrow(hits) > 0) {
          result_rows[[k]] <- data.frame(
            SNP = snp_id, GENE = hits$gene, VALUE = hits$value, STRATEGY = dir,
            stringsAsFactors = FALSE
          )
        } else {
          # Record NA if no match found
          result_rows[[k]] <- data.frame(
            SNP = snp_id, GENE = NA_character_, VALUE = NA_real_, STRATEGY = dir,
            stringsAsFactors = FALSE
          )
        }
      }
      # Combine results for the current chromosome
      do.call(rbind, result_rows)
    })

    # Combine results from all chromosomes for this strategy
    data_combination <- do.call(rbind, strategy_results_list)


    # --- 2. Preliminary Cleaning ---

    if (is.null(data_combination) || nrow(data_combination) == 0) {
      return(NULL)
    }

    colnames(data_combination) <- c("SNP", "GENE", "VALUE", "STRATEGY")
    data_combination$VALUE <- as.numeric(data_combination$VALUE)

    # Remove duplicates
    data_combination <- unique(data_combination)
    # Remove rows with no gene match
    data_combination <- data_combination[!is.na(data_combination$GENE), ]

    if (nrow(data_combination) == 0) {
      message(paste0("No valid links found for strategy: ", dir))
      return(NULL)
    }

    # --- 3. Gene Name Standardization & Filtering ---

    # Extract unique genes to minimize query time
    unique_genes <- unique(data_combination$GENE)

    # 3.1 Map Aliases to Official Symbols
    mapped_symbols <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique_genes,
      column = "SYMBOL",
      keytype = "ALIAS",
      multiVals = "first"
    ))

    # Retain original name if no mapping is found
    final_symbols <- ifelse(is.na(mapped_symbols), unique_genes, mapped_symbols)

    # 3.2 Retrieve Gene Types for Filtering
    gene_types <- suppressMessages(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = final_symbols,
      column = "GENETYPE",
      keytype = "SYMBOL",
      multiVals = "first"
    ))

    # 3.3 Create a Lookup Table
    gene_info_df <- data.frame(
      ORIGINAL_GENE = unique_genes,
      NEW_SYMBOL = final_symbols,
      GENE_TYPE = gene_types,
      stringsAsFactors = FALSE
    )

    # 3.4 Identify Protein-Coding Genes
    if (filter_gene) {
      # Strict Mode: Only keep protein-coding genes
      keep_mask <- !is.na(gene_info_df$GENE_TYPE) & gene_info_df$GENE_TYPE == "protein-coding"
    } else {
      # Permissive Mode: Keep all genes (just update their symbols)
      keep_mask <- rep(TRUE, nrow(gene_info_df))
    }

    # Subset valid mappings
    valid_map <- gene_info_df[keep_mask, ]

    # --- 4. Apply Mapping to Main Data ---

    # Inner join filters out invalid genes and updates names simultaneously
    df_cleaned <- merge(data_combination, valid_map, by.x = "GENE", by.y = "ORIGINAL_GENE")

    # Update GENE column with standardized symbols
    df_cleaned$GENE <- df_cleaned$NEW_SYMBOL
    df_cleaned <- df_cleaned[, c("SNP", "GENE", "VALUE", "STRATEGY")]

    # --- 5. Save Intermediate Results ---

    file_name <- paste0(pheno, "_", dir, "_raw_value")
    if(keep_intermediate_file) {
      save_path <- normalizePath(file.path(intermediate_file_folder, paste0(file_name, ".txt")))
      write.table(df_cleaned, file=save_path, sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
    }
    print(paste0(file_name, " successfully generated"))

    return(df_cleaned)
  })

  names(raw_results_list) <- dir_list
  raw_results_list <- raw_results_list[!sapply(raw_results_list, is.null)]


  #######################Step2: Intra-Strategy Normalization#######################################
  linking_score_list <- lapply(names(raw_results_list), function(strategy_name) {

    # 1. Retrieve the data for the current strategy
    df <- raw_results_list[[strategy_name]]

    # Check for empty dataframes
    if (is.null(df) || nrow(df) == 0) {
      return(NULL)
    }

    # 2. Apply Normalization Logic
    # Logic: For each SNP, keep only the gene(s) with the highest score.
    # If there are multiple genes with same max score, split the weight (1 / N).

    result_df <- df %>%
      dplyr::group_by(SNP) %>%
      dplyr::mutate(
        max_val = max(VALUE, na.rm = TRUE),
        is_max = (VALUE == max_val),
        count_max = sum(is_max, na.rm = TRUE)
      ) %>%
      dplyr::filter(is_max) %>%
      dplyr::mutate(VALUE = 1 / count_max) %>%
      dplyr::ungroup() %>%
      dplyr::select(SNP, GENE, VALUE, STRATEGY)

    # 3. Save Intermediate Results (Optional)
    if (keep_intermediate_file) {
      file_name <- paste0(pheno, "_", strategy_name, "_linking_score.txt")
      save_path <- normalizePath(file.path(intermediate_file_folder, file_name))

      write.table(result_df, file = save_path, sep = "\t",
                  row.names = FALSE, quote = FALSE, col.names = TRUE)
      print(paste0(file_name, " successfully generated"))
    }

    return(result_df)
  })

  # Restore the names of the list (lapply on names returns an unnamed list by default)
  names(linking_score_list) <- names(raw_results_list)

  # Remove any strategies that returned NULL (if any)
  linking_score_list <- linking_score_list[!sapply(linking_score_list, is.null)]

  # Validation check
  if (length(linking_score_list) == 0) {
    rlang::abort(paste("No valid linking scores generated for pheno:", pheno))
  }



  ##########################Step 3: combine the link result from seven strategy##########################

  # 1. Define Weights
  # Use default weights if custom_weights is not provided
  default_weights <- c(
    Exon = 100,
    Promoter = 100,
    finemappedciseQTLs_eQTLGen = 25,
    finemappedciseQTLs_GTeX = 7.5,
    EpiMap = 1.9,
    ABC = 1,
    Ciceroblood = 1
  )

  # Determine which weights to use
  use_weights <- if (is.numeric(custom_weights)) custom_weights else default_weights

  # 2. Merge Data
  # Combine the list from Step 2 into one dataframe
  if (length(linking_score_list) == 0) {
    rlang::abort(paste("No valid data to combine for pheno:", pheno))
  }
  combined_raw <- do.call(rbind, linking_score_list)

  # 3. Weighted Sum and Standardization

  final_df <- combined_raw %>%
    # Map weights to rows
    dplyr::mutate(weight_factor = use_weights[STRATEGY]) %>%

    # Calculate Weighted Sum for each SNP-Gene pair
    dplyr::group_by(SNP, GENE) %>%
    dplyr::summarise(
      RAW_SCORE = sum(VALUE * weight_factor, na.rm = TRUE),

      # Collapse strategy names
      STRATEGY_LIST = paste(unique(STRATEGY), collapse = ", "),

      .groups = 'drop'
    ) %>%

    # Clean strategy names
    dplyr::mutate(STRATEGY_LIST = gsub("ciseQTLs_", "", STRATEGY_LIST)) %>%

    # Calculate Fraction of Total SNP Score
    # Logic: For a given SNP, calculate the fraction of the total score belongs to this Gene
    dplyr::group_by(SNP) %>%
    dplyr::mutate(
      total_snp_score = sum(RAW_SCORE, na.rm = TRUE),

      # The Final Score: Ratio (Score / Total)
      WEIGHTED_SCORE = round(RAW_SCORE / total_snp_score, 3)
    ) %>%
    dplyr::ungroup()

  # 4. Export Result
  if (keep_intermediate_file) {
    save_path <- normalizePath(file.path(intermediate_file_folder, paste0(pheno, "_weighted_score.txt")))
    write.table(final_df, file = save_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
    print(paste0("Weighted score saved: ", save_path))
  }


  ##################Step 4: Precision Calculation##################

  # 1. Filter Low Confidence Links
  # Keep only links with weighted value > 0.5
  filtered_df <- final_df %>%
    dplyr::filter(WEIGHTED_SCORE > 0.5)

  if (nrow(filtered_df) == 0) {
    rlang::abort(paste("No links passed the threshold (> 0.5) for pheno:", pheno))
  }

  # 2. Calculate Combined Precision

  # Define precision dictionary
  default_precisions <- c(
    Exon = 1,
    finemappedGTeX = 0.676,
    EpiMap = 0.562,
    Promoter = 0.805,
    finemappedeQTLGen = 0.814,
    ABC = 0.469,
    Ciceroblood = 0.548
  )
  use_precisions <- if (is.numeric(custom_precisions)) custom_precisions else default_precisions

  # calculate precision: P = 1 - product(1 - p_i)
  calc_precision_row <- function(strategy_str, precision_map) {
    # Split string by comma (e.g., "Exon, Promoter")
    strat_vec <- unlist(strsplit(strategy_str, ",\\s*"))

    # Map strategies to precision values
    p_vals <- precision_map[strat_vec]

    # Handle missing/NA values (if strategy name doesn't match dict)
    if (any(is.na(p_vals))) {
      warning(paste("Missing precision for strategy:", paste(strat_vec[is.na(p_vals)], collapse=", ")))
      p_vals <- na.omit(p_vals)
    }

    if (length(p_vals) == 0) return(0)

    combined_val <- 1 - prod(1 - p_vals)
    return(round(combined_val, 3))
  }

  # Apply precision calculation
  filtered_df$PRECISION <- sapply(filtered_df$STRATEGY_LIST, calc_precision_row, precision_map = use_precisions)

  # 3. Calculate Confidence Score (Combine with PIP)

  if (causal_SNP) {
    # Ensure SNP_data has 'SNP' and 'PIP' columns
    pip_info <- SNP_data %>% dplyr::select(SNP, PIP)

    filtered_df <- filtered_df %>%
      dplyr::left_join(pip_info, by = "SNP") %>%
      dplyr::mutate(CONFIDENCE_SCORE = PRECISION * PIP)
  } else {
    # If no causal SNP data, use assumed PIP
    current_pip <- if (!is.null(assumed_PIP)) assumed_PIP else 0.5
    filtered_df <- filtered_df %>%
      dplyr::mutate(CONFIDENCE_SCORE = PRECISION * current_pip)
  }


  # 4. Map Symbol to ENSGID
  genes_to_map <- unique(filtered_df$GENE)

  ens_ids <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = genes_to_map,
    column = "ENSEMBL",
    keytype = "SYMBOL",
    multiVals = "first"
  ))

  total_genes <- length(genes_to_map)
  mapped_genes <- sum(!is.na(ens_ids))
  mapping_rate <- round((mapped_genes / total_genes) * 100, 2)

  rlang::inform(
    message = c(
      paste0("Gene annotation completed (Symbol -> ENSGID)."),
      "i" = paste0("Mapping success rate: ", mapping_rate, "% (", mapped_genes, "/", total_genes, " genes)."),
      if (mapping_rate < 100) "i" = "Unmapped genes are set to NA." else NULL
    )
  )

  annot_lookup <- data.frame(
    GENE = genes_to_map,
    ENSGID = ens_ids,
    stringsAsFactors = FALSE
  )

  filtered_df <- filtered_df %>%
    dplyr::left_join(annot_lookup, by = "GENE") %>% dplyr::select(-total_snp_score) %>%
    dplyr::rename(SYMBOL = GENE,
                  STRATEGY = STRATEGY_LIST) %>%
    dplyr::select(SNP,SYMBOL,ENSGID,RAW_SCORE, WEIGHTED_SCORE, STRATEGY, PRECISION, PIP, CONFIDENCE_SCORE)



  # --- 5. Export Final Result ---
  out_file <- normalizePath(file.path(save_folder, paste0(pheno, "_cS2G.csv")))
  write.csv(filtered_df, out_file, row.names = FALSE)

  print(paste0("Process completed. Result saved to: ", out_file))

  return(filtered_df)
}
