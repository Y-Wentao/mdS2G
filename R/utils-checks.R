#' Validate Phenotype Name
#' @noRd
validate_pheno_name <- function(pheno_name) {
  if (missing(pheno_name) || !is.character(pheno_name) || length(pheno_name) == 0) {
    rlang::abort(
      message = c(
        'Invalid "pheno_name" input.',
        "i" = 'Please provide a character vector, e.g., c("height", "BMI").'
      )
    )
  }

  pheno_name <- trimws(pheno_name)
  pheno_name <- pheno_name[pheno_name != ""]
  pheno_name <- unique(pheno_name)
  if (length(pheno_name) == 0) {
    rlang::abort("'pheno_name' contains no valid values.")
  }
  return(pheno_name)
}


#' Process and Validate SNP Data
#' @importFrom data.table fread
#' @noRd
process_snp_data <- function(SNP_data, causal_SNP, pheno_name) {
  if (missing(SNP_data)) {
    rlang::abort("'SNP_data' must be provided.")
  }

  # 1. Load data
  if (is.character(SNP_data) && length(SNP_data) == 1) {
    if (!file.exists(SNP_data)) rlang::abort(paste("SNP_data file not found:", SNP_data))
    SNP_data <- as.data.frame(data.table::fread(SNP_data, header = TRUE, data.table = FALSE))
  }

  if (!is.data.frame(SNP_data)) {
    rlang::abort("'SNP_data' must be a data frame or a valid file path.")
  }

  # 2. Standardize column names to Upper Case
  # This allows "chr", "Chr", "CHR" but rejects "chrom", "pos"
  colnames(SNP_data) <- toupper(colnames(SNP_data))

  # 3. Define Required Columns
  required_cols <- c("CHR", "SNP", "BP")

  if (causal_SNP) {
    required_cols <- c(required_cols, "PIP")
  }

  # Only require PHENO column if multiple phenotypes are requested
  if (length(pheno_name) > 1) {
    required_cols <- c(required_cols, "PHENO")
  }

  # 4. Check Missing Columns
  missing_cols <- setdiff(required_cols, colnames(SNP_data))
  if (length(missing_cols) > 0) {
    rlang::abort(
      c(
        "Missing required columns in 'SNP_data'.",
        "x" = paste0("Missing columns: ", paste(missing_cols, collapse = ", ")),
        "i" = "Columns must be explicitly named: CHR, SNP, BP (and PIP, PHENO if applicable)."
      )
    )
  }

  # 5. Subset to keep only relevant columns (and ensure correct order)
  if (!"PHENO" %in% colnames(SNP_data)) {
    SNP_data$PHENO <- pheno_name[1]
  }

  # 6. Validate
  final_cols <- c("CHR", "SNP", "BP", "PHENO")
  if (causal_SNP) final_cols <- c(final_cols, "PIP")
  SNP_data <- SNP_data[, final_cols, drop = FALSE]
  unique_pheno_in_data <- unique(SNP_data$PHENO)
  if (!all(unique(pheno_name) %in% unique_pheno_in_data)) {
    unmatched <- setdiff(pheno_name, unique_pheno_in_data)
    rlang::abort(paste("The following 'pheno_name' entries not found in SNP_data:", paste(unmatched, collapse = ", ")))
  }

  # 7. Standardize CHR values (Data Cleaning)
  # Even with correct column names, content might still have "chr1" instead of "1"
  SNP_data$CHR <- gsub("chr", "", SNP_data$CHR, ignore.case = TRUE)
  SNP_data$CHR <- as.numeric(trimws(SNP_data$CHR))

  # 8. Filter Autosomal & Validate Row Counts (NEW LOGIC)
  SNP_data <- SNP_data[SNP_data$CHR %in% 1:22, ]

  if (nrow(SNP_data) == 0) {
    rlang::abort(
      c(
        "No valid autosomal SNPs remained after filtering.",
        "x" = "All rows were removed because CHR values were not in range 1-22.",
        "i" = "Please check the 'CHR' column content. Only autosomal chromosomes (1-22) are supported."
      )
    )
  } else {
    # Inform user about the valid input count
    message(sprintf("Total valid autosomal variants input: %d. Please check if this meets expectations.", nrow(SNP_data)))
  }

  return(SNP_data)
}


#' Validate S2G Folder Structure
#' @noRd
validate_s2g_folder <- function(folder) {
  if (missing(folder) || !is.character(folder) || length(folder) != 1 || !dir.exists(folder)) {
    rlang::abort("'s2g_folder' must be provided and must be a valid directory path.")
  }
  folder <- normalizePath(folder)
  required_subfolders <- c("Exon", "Promoter", "finemappedciseQTLs_eQTLGen",
                           "finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood")

  missing_subs <- required_subfolders[!dir.exists(file.path(folder, required_subfolders))]
  if (length(missing_subs) > 0) {
    rlang::abort(paste("Missing s2g subfolders:", paste(missing_subs, collapse = ", ")))
  }

  # Check file existence inside subfolders
  for (sub in required_subfolders) {
    for (chr in 1:22) {
      bed_file <- file.path(folder, sub, paste0("chr", chr, ".bed.gz"))
      if (!file.exists(bed_file)) {
        rlang::abort(paste0("Missing annotation file: ", bed_file))
      }
    }
  }
  return(folder)
}

#' Validate Writeable Directory
#' @noRd
validate_dir_writeable <- function(path, arg_name) {
  if (missing(path) || !is.character(path) || length(path) != 1 || !dir.exists(path)) {
    rlang::abort(paste0("'", arg_name, "' must be a valid directory path."))
  }
  if (file.access(path, 2) != 0) {
    rlang::abort(paste0("User does not have write permissions for '", arg_name, "'."))
  }
  return(normalizePath(path))
}

#' Validate Assumed PIP
#' @noRd
validate_assumed_pip <- function(causal_SNP, assumed_PIP) {
  if (!causal_SNP) {
    if (is.null(assumed_PIP)) {
      rlang::warn("assumed_PIP not provided, defaulting to 0.5.")
      return(0.5)
    }
    if (!is.numeric(assumed_PIP)) rlang::abort("`assumed_PIP` must be numeric.")
    return(assumed_PIP)
  } else {
    if (!is.null(assumed_PIP)) {
      rlang::abort("`causal_SNP` must be FALSE when `assumed_PIP` is provided.")
    }
    return(NULL)
  }
}


#' Process Weights or Precisions
#' @param input User input (FALSE, NULL, or file path).
#' @param arg_name Name of the argument for error reporting.
#' @return A named numeric vector or NULL.
#' @importFrom data.table fread
#' @importFrom rlang abort
#' @noRd
process_weights_or_precisions <- function(input, arg_name) {

  # 1. Check for default usage (NULL or FALSE)
  if (is.null(input) || identical(input, FALSE)) {
    return(NULL)
  }

  # 2. Strict Type Check: Input must be a single character string (file path)
  if (!is.character(input) || length(input) != 1) {
    rlang::abort(
      message = c(
        paste0("Invalid input for '", arg_name, "'."),
        "i" = "It must be either FALSE (to use defaults) or a valid file path string.",
        "x" = paste0("You provided an object of type: ", class(input)[1])
      )
    )
  }

  # 3. Verify file existence
  if (!file.exists(input)) {
    rlang::abort(paste0("File not found for '", arg_name, "': ", input))
  }

  # 4. Attempt to read the file
  data <- tryCatch({
    data.table::fread(input, header = FALSE, data.table = FALSE)
  }, error = function(e) {
    rlang::abort(paste0("Failed to read file '", input, "'. Ensure it is a valid CSV/TXT."))
  })

  # 5. Validate dimensions (Minimum 7 strategies, 2 columns)
  if (nrow(data) < 7 || ncol(data) < 2) {
    rlang::abort(paste0("The file for '", arg_name, "' must have at least 7 rows and 2 columns."))
  }

  # 6. Convert to named numeric vector (Col 1 = Names, Col 2 = Values)
  vec <- setNames(as.numeric(data[[2]]), data[[1]])

  # 7. Validate required keys
  required_keys <- c("Exon", "Promoter", "finemappedciseQTLs_eQTLGen",
                     "finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood")

  missing_keys <- setdiff(required_keys, names(vec))
  if (length(missing_keys) > 0) {
    rlang::abort(
      message = c(
        paste0("Missing required strategies in '", arg_name, "' file."),
        "x" = paste0("Missing: ", paste(missing_keys, collapse = ", "))
      )
    )
  }

  # 8. Validate numeric values
  if (any(is.na(vec))) {
    rlang::abort(paste0("The second column of '", arg_name, "' must contain valid numeric values."))
  }

  return(as.numeric(vec))
}

#' Check Directories
#' @noRd
check_dir <- function(path, name) {
  if (missing(path) || !is.character(path) || length(path) != 1 || !dir.exists(path)) {
    rlang::abort(paste0("'", name, "' must be provided and must be a valid directory path."))
  }
  return(normalizePath(path, mustWork = TRUE))
}
