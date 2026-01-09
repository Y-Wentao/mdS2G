#' @title Link Causal SNPs to Genes using cS2G Strategy
#'
#' @description
#' This function implements the **locus-based linking** step of the \strong{mdS2G framework}.
#' It adapts the published combined SNP-to-Gene (cS2G) strategy to link fine-mapped causal SNPs
#' to their putative target genes by integrating evidence from seven sub-strategies.
#' A \strong{Confidence Score} is assigned to each SNP-Gene pair by weighting the fine-mapping
#' posterior probability with the linking precision.
#'
#' Since the cS2G annotation resources rely on physical position-based mapping, and the underlying annotation files are built on a specific assembly,
#' the physical positions in the input `SNP_data` \strong{MUST be based on the hg19 (GRCh37) reference genome}.
#' Using coordinates from other builds (e.g., hg38) will result in incorrect gene assignments.
#'
#' @details
#' The function operates based on the cS2G framework but introduces two key modifications
#' tailored for the mdS2G workflow:
#'
#' \itemize{
#'   \item \strong{Multiplicative Precision Update:} For SNP-Gene pairs supported by
#'   multiple sub-strategies simultaneously, we employ a \strong{Bayesian integration approach}
#'   to calculate a joint posterior probability. Assuming independence among the different
#'   functional genomic evidence layers, this method combines individual precisions to
#'   reinforce the confidence score when convergent evidence is present.
#'
#'   \item \strong{Inclusive Retention Strategy:} Links with a Confidence Score \eqn{\ge 0.8}
#'   are flagged as \strong{credible}. Unlike standard thresholds that discard data,
#'   this function \strong{retains all identified links} (including low-confidence ones)
#'   to allow for potential validation by subsequent evidence layers in the framework.
#' }
#'
#' @param pheno_name The pheno or trait include in your study, for example "height"; If you have two or more pheno, please separate them with commas, for example "height,BMI".
#' @param SNP_data The path point to your causal SNP or lead SNP (it is not recommended yet) result, but you should never use them together. The data should have each line representing a SNP, and it should contain at least CHR, SNP, BP, PHENO (if multiple phenotypes are provided in pheno_name), and PIPS (if they are causal SNPs).
#' @param s2g_folder The path point to 'S2G_bed' mapping file. This file should be compressed and contain "Exon", "Promoter", "finemappedciseQTLs_eQTLGen", finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood" seven subfolders. The specific files within these subfolders should be uncompressed
#' @param save_folder The path where you want to store the result file.
#' @param keep_intermediate_file Whether you want to store intermediate files, default to False.
#' @param intermediate_file_folder If you want to store intermediate files, you should provide a specific storage path. It's default to NULL
#' @param filter_gene We suggest deleting non coding genes to clarify the results. If you do not want to exclude, please select False. It's default to TRUE.
#' @param custom_weights We used the optimal weights provided in the Gazal manuscript as the integrated weights for the seven strategies. If you have any other weights you would like to use, please provide the specific file. This data file should not have a header. The first column corresponds to the names of the seven strategies used (Please refer to the description of the s2d_folder parameter above for details), the second column should correspond to specific weights. That is to say, this file should have seven rows and two columns. It's default to NULL.
#' @param custom_precisions We used the precision of these seven strategies from the Gazal manuscript in the validation set as part of the confidence score calculation. If you have any other data you would like to use, please enter it. The specific requirements of the document are the same as' other weights'. It's default to NULL.
#' @param causal_SNP Whether the SNPs included in 'SNP_data' are causal SNPs, default to TRUE.
#' @param assumed_PIP If your SNP does not come from fine mapping, you need to provide a hypothetical posterior probability to calculate the confidence score. We recommend using a conservative data (default is 0.5), although this may require support from the other two types of evidence for all linking results, it can effectively reduce false positive results.
#' @param keep_unmatched_gene Whether to retain genes that failed to convert ENSGID, default to False
#'
#' @return A data.frame and a file named pheno_cS2G.csv
#' @export

cs2g_link <- function(pheno_name,
                      SNP_data,
                      s2g_folder,
                      save_folder,
                      keep_intermediate_file = FALSE,
                      intermediate_file_folder = NULL,
                      filter_gene = TRUE,
                      custom_weights = FALSE,
                      custom_precisions = FALSE,
                      causal_SNP = TRUE,
                      assumed_PIP = NULL,
                      keep_unmatched_gene = FALSE) {

  # 1. Standardize Phenotype Name
  pheno_name <- validate_pheno_name(pheno_name)
  rlang::inform(
    message = c(
      "i" = paste0("Phenotypes to be analyzed (", length(pheno_name), "):"),
      " " = paste(pheno_name, collapse = ", ")
    )
  )

  # 2. Validate Flags
  if (!is.logical(causal_SNP) || length(causal_SNP) != 1) {
    rlang::abort("'causal_SNP' must be TRUE or FALSE.")
  }
  if (!is.logical(keep_intermediate_file) || length(keep_intermediate_file) != 1) {
    rlang::abort("'keep_intermediate_file' must be TRUE or FALSE.")
  }
  if (!is.logical(keep_unmatched_gene) || length(keep_unmatched_gene) != 1) {
    rlang::abort("'keep_unmatched_gene' must be TRUE or FALSE.")
  }
  if (!is.logical(filter_gene) || length(filter_gene) != 1) {
    rlang::abort("'filter_gene' must be TRUE or FALSE.")
  }


  # 3. Assumed PIP needs to be provided if not a fine-mapping result
  if (!causal_SNP) {
    if (is.null(assumed_PIP)) {
      rlang::warn("assumed_PIP not provided, defaulting to 0.5.")
      assumed_PIP <- 0.5
    }
    if (!is.numeric(assumed_PIP)) rlang::abort("`assumed_PIP` must be numeric.")
  } else {
    if (!is.null(assumed_PIP)) {
      rlang::abort("`causal_SNP` must be FALSE when `assumed_PIP` is provided.")
    }
  }

  # 4. Process SNP Data
  #    Returns a standardized data frame
  SNP_data <- process_snp_data(SNP_data, causal_SNP, pheno_name)


  # 5. Validate Paths of Folders
  s2g_folder <- validate_s2g_folder(s2g_folder)
  save_folder <- validate_dir_writeable(save_folder, "save_folder")


  # 6. Handle Intermediate Files logic
  if (keep_intermediate_file) {
    if (is.null(intermediate_file_folder)) {
      rlang::abort("'intermediate_file_folder' must be provided when 'keep_intermediate_file' is TRUE.")
    }
    intermediate_file_folder <- validate_dir_writeable(intermediate_file_folder, "intermediate_file_folder")
  } else if (!is.null(intermediate_file_folder)) {
    rlang::abort("'keep_intermediate_file' must be TRUE when 'intermediate_file_folder' is provided.")
  }

  # 7. Process Auxiliary Data Files
  custom_weights_vec    <- process_weights_or_precisions(custom_weights, "custom_weights")
  custom_precisions_vec <- process_weights_or_precisions(custom_precisions, "custom_precisions")

  # --- End of Validation ---
  message("Inputs validated successfully. Starting analysis...")


  ##############################cS2G-mapping#################################
  for (pheno in pheno_name) {
    message(paste0("\n>>> Processing Phenotype: ", pheno))

    # Data Slicing
    current_SNP_data <- SNP_data[SNP_data$PHENO == pheno, ]

    if (nrow(current_SNP_data) == 0) {
      warning(paste("No SNP rows found for phenotype:", pheno, "- Skipping."))
      next
    }

    tryCatch({
      link_cs2g(
        pheno_name = pheno,
        SNP_data = current_SNP_data,
        s2g_folder = s2g_folder,
        save_folder = save_folder,
        keep_intermediate_file = keep_intermediate_file,
        intermediate_file_folder = intermediate_file_folder,
        filter_gene = filter_gene,
        custom_weights = custom_weights_vec,
        custom_precisions = custom_precisions_vec,
        causal_SNP = causal_SNP,
        assumed_PIP = assumed_PIP,
        keep_unmatched_gene = keep_unmatched_gene
      )

      message(paste(">>> Successfully finished:", pheno))

    }, error = function(e) {
      message(paste("!!! Error processing phenotype", pheno, ":"))
      message(e$message)
      warning(paste("Skipping", pheno, "due to calculation error."))
    })
  }
}
