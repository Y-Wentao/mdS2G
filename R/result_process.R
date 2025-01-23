#' @title Process results from mdS2G
#'
#' @description This function is designed to process and summarize the results previously obtained from mdS2G.
#'
#' @param coloc_result
#' @param SNP_data_create_above
#' @param cs2g_cred_data
#' @param save_folder
#'
#' @return
#' @export
#'
#' @examples

result_process <- function(coloc_result, SNP_data_create_above, cs2g_cred_data, save_folder)
{
  suppressMessages(suppressWarnings(library(dplyr)))
  suppressMessages(suppressWarnings(library(stringr)))

  if (missing(coloc_result) || !is.character(coloc_result) || length(coloc_result) != 1 || !file.exists(coloc_result)) {
    rlang::abort("'coloc_result' must be provided as a file path.")
  }
  coloc_result <- data.table::fread(coloc_result, stringsAsFactors = FALSE, data.table=FALSE)
  if (!is.data.frame(coloc_result) && !is.matrix(coloc_result)) {
    rlang::abort("'coloc_result' must be a valid file path or a data frame/matrix.")
  }
  expected_colnames <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "PHENO", "SNP", "GENE", "TISSUE")
  actual_colnames <- colnames(coloc_result)
  if (!identical(actual_colnames, expected_colnames)) {
    rlang::abort(paste0(
      "The expected coloc_result file column names should be:\n",
      paste(expected_colnames, collapse = "\t"), "\n",
      "But the provided file column names are:\n",
      paste(actual_colnames, collapse = "\t"), "\n",
      "Please check the input file for any issues."
    ))
  }

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



  if (missing(cs2g_cred_data) || !is.character(cs2g_cred_data) || length(cs2g_cred_data) != 1 || !file.exists(cs2g_cred_data)) {
    rlang::abort("'cs2g_cred_data' must be provided as a file path.")
  }
  cs2g_cred_data <- data.table::fread(cs2g_cred_data, stringsAsFactors = FALSE, data.table=FALSE)
  if (!is.data.frame(cs2g_cred_data) && !is.matrix(cs2g_cred_data)) {
    rlang::abort("'cs2g_cred_data' must be a valid file path or a data frame/matrix.")
  }
  expected_colnames <- c("SNP", "GENE", "VALUE", "STRATEGY", "VALUE_WEIGHTED", "precision", "PIP", "confidence_score", "ENSGID", "PHENO")
  actual_colnames <- colnames(cs2g_cred_data)
  if (!identical(actual_colnames, expected_colnames)) {
    rlang::abort(paste0(
      "The expected cs2g_cred_data file column names should be:\n",
      paste(expected_colnames, collapse = "\t"), "\n",
      "But the provided file column names are:\n",
      paste(actual_colnames, collapse = "\t"), "\n",
      "Please check the input file for any issues."
    ))
  }

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



  ##########merge data
  coloc_result$PHENOSNP <- paste0(coloc_result$PHENO,"-",coloc_result$SNP)
  data_wide <- coloc_result %>%
    dplyr::group_by(PHENOSNP, GENE) %>%
    dplyr::summarise(TISSUE = paste(TISSUE, collapse = ","), .groups = 'drop') %>%
    dplyr::group_by(PHENOSNP) %>%
    dplyr::mutate(row = row_number()) %>%
    tidyr::pivot_wider(
      names_from = row,
      values_from = c(GENE, TISSUE),
      names_sep = ""
    ) %>%
    ungroup() %>%
    tidyr::separate(PHENOSNP, into = c("PHENO", "SNP"), sep = "-") %>%
    rename_with(~sub("TISSUE", "TISSUE_GENE", .), starts_with("TISSUE"))

  data_long <- data_wide %>%
    tidyr::pivot_longer(cols = starts_with("GENE"), names_to = "GENE_col", values_to = "GENE") %>%
    tidyr::pivot_longer(cols = starts_with("TISSUE_GENE"), names_to = "TISSUE_col", values_to = "TISSUE") %>%
    dplyr::filter(str_extract(GENE_col, "\\d+") == str_extract(TISSUE_col, "\\d+")) %>%
    dplyr::filter(!is.na(GENE) & !is.na(TISSUE))


  unique_pairs <- data_long %>%
    distinct(PHENO, SNP, GENE) %>%
    nrow()


  avg_tissue_per_pair <- data_long %>%
    group_by(PHENO, SNP, GENE) %>%
    summarise(average_tissue_count = mean(str_count(TISSUE, ",") + 1), .groups = "drop")


  message(paste0("Total number of SNP-GENE-PHENO pairs with colocalization evidence: ", unique_pairs, "\n"))
  message(paste0("Average number of tissues supporting each colocalization pair: ", mean(avg_tissue_per_pair$average_tissue_count)))
  output_content <- paste0(
    "Total number of SNP-GENE-PHENO pairs with colocalization evidence: ", unique_pairs, "\n",
    "Average number of tissues supporting each colocalization pair: ", mean(avg_tissue_per_pair$average_tissue_count), "\n"
  )
  writeLines(output_content, file.path(save_folder,"coloc_metadata.txt"))
  write.csv(data_wide,file.path(save_folder,"coloc_result_processed.csv"),row.names = F)



  colnames(SNP_data) <- colnames(SNP_data) %>%
    stringr::str_replace_all("^GENE([1-5])$", "POPS_GENE_\\1") %>%
    stringr::str_replace("^GENE6$", "cS2G_GENE")
  SNP_data$PHENOSNP <- paste0(SNP_data$PHENO,"-",SNP_data$SNP)

  coloc_data <- data_wide
  coloc_data$PHENOSNP <- paste0(coloc_data$PHENO,"-",coloc_data$SNP)
  colnames(coloc_data) <- colnames(coloc_data) %>%
    stringr::str_replace_all("^GENE([1-5])$", "COLOC_GENE_\\1") %>%
    stringr::str_replace_all("^TISSUE_GENE([1-5])$", "COLOC_TISSUE_\\1")

  coloc_columns <- coloc_data %>%
    select(starts_with("COLOC_GENE"), starts_with("COLOC_TISSUE")) %>%
    colnames()

  merged_data <- SNP_data %>%
    left_join(coloc_data[, c("PHENOSNP", coloc_columns)], by = "PHENOSNP") %>%
    select(-PHENOSNP) %>% dplyr::mutate_all(~tidyr::replace_na(., "NONE"))


  base_data <- merged_data %>% select(PHENO, CHR, SNP, BP)

  coloc_gene_columns <- grep("^COLOC_GENE_", colnames(merged_data), value = TRUE)
  for (i in 1:nrow(merged_data)) {
    genes <- unlist(merged_data[i, c("POPS_GENE_1", "POPS_GENE_2", "POPS_GENE_3", "POPS_GENE_4", "POPS_GENE_5", "cS2G_GENE", coloc_gene_columns)], use.names = FALSE)

    genes <- genes[genes != "NONE"]

    if (length(genes) == 0) {
      next
    }

    gene_evidence <- rep("", length(genes))
    names(gene_evidence) <- genes

    for (j in 1:5) {
      gene <- merged_data[i, paste0("POPS_GENE_", j)]
      if (gene != "NONE") {
        gene_evidence[gene] <- paste(gene_evidence[gene], "POPS", sep = ifelse(gene_evidence[gene] == "", "", ","))
      }
    }

    if (merged_data[i, "cS2G_GENE"] != "NONE") {
      gene <- merged_data[i, "cS2G_GENE"]
      gene_evidence[gene] <- paste(gene_evidence[gene], "cS2G", sep = ifelse(gene_evidence[gene] == "", "", ","))
    }

    for (j in 1:length(coloc_gene_columns)) {
      gene <- merged_data[i, paste0("COLOC_GENE_", j)]
      if (gene != "NONE") {
        gene_evidence[gene] <- paste(gene_evidence[gene], "coloc", sep = ifelse(gene_evidence[gene] == "", "", ","))
      }
    }

    gene_count <- table(genes)

    supported_genes <- names(gene_count[gene_count >= 2])

    if (length(supported_genes) > 0) {
      for (j in 1:length(supported_genes)) {
        base_data[i, paste0("CREDIBLE_GENE_", j)] <- supported_genes[j]
        base_data[i, paste0("EVIDENCE_GENE_", j)] <- gene_evidence[supported_genes[j]]
      }
    }
  }

  data_muti_envi <- base_data
  write.csv(data_muti_envi,file.path(save_folder,"credible_gene_muti_evidence.csv"),row.names=F)


  cs2g_cred_data$EVIDENCE_GENE_highConfi <- "cS2G(confi>=0.8)"
  colnames(cs2g_cred_data)[colnames(cs2g_cred_data)=="ENSGID"] <- "CREDIBLE_GENE_highConfi"


  data_muti_envi$PHENOSNP <- paste0(data_muti_envi$PHENO,"-",data_muti_envi$SNP)
  cs2g_cred_data$PHENOSNP <- paste0(cs2g_cred_data$PHENO,"-",cs2g_cred_data$SNP)


  merged_data <- data_muti_envi %>%
    left_join(cs2g_cred_data[,c("PHENOSNP","CREDIBLE_GENE_highConfi","EVIDENCE_GENE_highConfi")],by="PHENOSNP") %>% select(-PHENOSNP)

  write.csv(merged_data, file.path(save_folder,"CredibleGene_all.csv"), row.names=F)


  merged_clean <- merged_data %>%
    tidyr::pivot_longer(cols = starts_with("CREDIBLE_GENE_"), names_to = "Gene_Type", values_to = "Gene") %>%
    dplyr::filter(!is.na(Gene)) %>%
    dplyr::select(PHENO, Gene)


  phenotype_genes <- list()

  phenotypes <- unique(merged_data$PHENO)
  for (pheno in phenotypes) {
    unique_genes <- merged_clean %>%
      dplyr::filter(PHENO == pheno) %>%
      dplyr::pull(Gene) %>%
      unique()

    phenotype_genes[[pheno]] <- unique_genes
  }


  max_length <- max(sapply(phenotype_genes, length))


  unique_gene <- as.data.frame(matrix(NA, nrow = max_length, ncol = length(phenotypes)))
  colnames(unique_gene) <- phenotypes


  for (pheno in phenotypes) {
    unique_gene[[pheno]][1:length(phenotype_genes[[pheno]])] <- phenotype_genes[[pheno]]
  }
  write.csv(unique_gene, file.path(save_folder,"unique_CredibleGene.csv"),row.names=F)

}







