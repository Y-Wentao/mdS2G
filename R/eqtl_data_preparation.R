#' Prepare eQTL summary statistics for colocalization
#'
#' @description
#' Due to its high accessibility and standardized formatting, this function is specifically designed to process GTEx v8 (accessible at \url{https://yanglab.westlake.edu.cn/software/smr/#DataResource}) and eQTLGEN (accessible at in \url{https://www.eqtlgen.org/cis-eqtls.html}) cis-eQTL summary statistics in **SMR binary format**.
#'
#' \strong{Important Notes:}
#' \itemize{
#'   \item \strong{Data Completeness:} Do \strong{not} use "Lite" summary statistics (which only contain significant SNPs reaching a certain threshold). Colocalization analysis requires \strong{full summary statistics} (including non-significant SNPs) to accurately evaluate the posterior probability of a shared causal variant.
#'   \item \strong{Other Data Sources:} While optimized for YangLab's format, you may use GTEx, eQTLGen, or other QTL summary statistics from different sources. However, if bypassing this function's processing, you must manually format your data into individual text files per gene, matching the output structure of this function. Matching relies on names here, so the downloaded summary file should be in compressed format and the name should not be changed.
#'   \item \strong{MAF Imputation:} The SMR binary files provided by YangLab lack Minor Allele Frequency (MAF) information. This function automatically imputes MAF using the **GTEx V8 reference data** bundled with this package. (Note: eQTLGen data does not require imputation).
#' }
#'
#' @param eqtl_data_folder Directory containing the raw eQTL data.
#' @param eqtl_data_tissue Character vector. A vector of specific tissue names to process (e.g., \code{c("Whole_Blood", "Liver")}). If \code{NULL} (default), the function detects and processes all available tissues in the folder. Note: Input names are validated against the standard GTEx V8/eQTLGen tissue list.
#' @param save_folder Directory to save output. Defaults to eqtl_data_folder.
#' @param smr_exe Path to the SMR executable file.
#' @param extract_gene_list Path to a file containing a list of genes to extract.
#' @param smr_file_delete Logical. Whether to delete intermediate SMR binary files after conversion.
#' @param complete_file_delete Logical. Whether to delete the huge converted text file after splitting.
#' @param maf_impute Logical. Whether to impute MAF from an external reference.
#' @param maf_impute_file Path to the external MAF reference file (containing 'SNP' and 'MAF' columns). By default, the function uses the GTEx V8 MAF reference data provided internally by this package.
#'
#' @export
eqtl_data_preparation <- function(eqtl_data_folder,
                                  eqtl_data_tissue = NULL,
                                  save_folder = NULL,
                                  smr_exe,
                                  extract_gene_list,
                                  smr_file_delete = TRUE,
                                  complete_file_delete = TRUE,
                                  maf_impute = TRUE,
                                  maf_impute_file = NULL) {

  # --- 1. Dependency & Environment Setup ---
  # Valid tissue list
  valid_tissue_names <- c(
    "Thyroid", "Adipose_Visceral_Omentum", "Artery_Aorta", "Kidney_Cortex", "Heart_Left_Ventricle",
    "Small_Intestine_Terminal_Ileum", "Skin_Sun_Exposed_Lower_leg", "Pancreas", "Pituitary",
    "Brain_Caudate_basal_ganglia", "Minor_Salivary_Gland", "Brain_Cerebellar_Hemisphere",
    "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Vagina", "Artery_Tibial", "Muscle_Skeletal",
    "Brain_Anterior_cingulate_cortex_BA24", "Ovary", "Esophagus_Muscularis", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Substantia_nigra", "Testis",
    "Adipose_Subcutaneous", "Artery_Coronary", "Adrenal_Gland", "Cells_Cultured_fibroblasts",
    "Cells_EBV-transformed_lymphocytes", "Lung", "Uterus", "Brain_Amygdala", "Nerve_Tibial", "Spleen",
    "Stomach", "Colon_Sigmoid", "Esophagus_Mucosa", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction",
    "Heart_Atrial_Appendage", "Brain_Cerebellum", "Brain_Cortex", "Brain_Hypothalamus",
    "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Breast_Mammary_Tissue", "Liver",
    "Whole_Blood", "eQTLGEN"
  )

  # --- 2. Input Validation ---
  # Check Folders
  eqtl_data_folder <- check_dir(eqtl_data_folder, "eqtl_data_folder")
  eqtl_data_folder <- normalizePath(eqtl_data_folder)

  if (is.null(save_folder)) {
    save_folder <- eqtl_data_folder
  } else {
    if (!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
    if (file.access(save_folder, 2) != 0) rlang::abort("'save_folder' is not writable.")
    save_folder <- normalizePath(save_folder)
  }

  # Check SMR Executable (Strict Requirement)
  if (missing(smr_exe) || !file.exists(smr_exe)) {
    rlang::abort("Error: 'smr_exe' path is invalid. SMR tool is required for this function.")
  }

  # Check Gene List
  if (missing(extract_gene_list) || !file.exists(extract_gene_list)) rlang::abort("'extract_gene_list' not found.")
  target_genes <- as.character(data.table::fread(extract_gene_list, header = FALSE)[[1]])

  # Load MAF Data
  maf_data <- NULL
  if (maf_impute) {
    if (is.null(maf_impute_file)) maf_impute_file <- system.file("extdata", "maf_impute.csv.gz", package = "mds2g")
    if (!file.exists(maf_impute_file)) rlang::abort("'maf_impute_file' not found.")
    maf_data <- data.table::fread(maf_impute_file, select = c("SNP", "MAF"), data.table = FALSE)
  }

  # --- 2. Tissue Detection (Strict SMR Zip Mode) ---
  zip_pattern <- "\\.zip$|cis-eQTL-SMR_20191212\\.tar\\.gz$"
  detected_files <- list.files(eqtl_data_folder, pattern = zip_pattern, full.names = FALSE)

  if (length(detected_files) == 0) {
    rlang::abort(paste0("No valid compressed SMR files (.zip or .tar.gz) found in: ", eqtl_data_folder))
  }

  # Map filenames to Tissue IDs
  detected_tissues <- gsub("\\.zip$|\\.tar\\.gz$", "", detected_files)
  detected_tissues[detected_tissues == "cis-eQTL-SMR_20191212"] <- "eQTLGEN"

  # Validate User Input vs Detected Files
  if (is.null(eqtl_data_tissue)) {
    tissue_list <- detected_tissues
    message("Detected ", length(tissue_list), " compressed files for processing.")
  } else {
    # Standardization
    eqtl_data_tissue <- trimws(unlist(strsplit(eqtl_data_tissue, ",")))
    missing_files <- setdiff(eqtl_data_tissue, detected_tissues)
    if (length(missing_files) > 0) {
      rlang::abort(paste0("Requested tissues not found ", paste(missing_files, collapse = ", ")))
    }
    tissue_list <- eqtl_data_tissue
  }

  # Warn about non-standard names (Informational only)
  invalid_tissues <- setdiff(tissue_list, valid_tissue_names)
  if (length(invalid_tissues) > 0) {
    warning("Some tissues are not in the standard GTEx/eQTLGen list: ", paste(invalid_tissues, collapse = ", "), call. = FALSE)
  }


  # --- 3. Helper Functions ---

  # A. SMR Conversion
  run_smr <- function(bfile, outfile) {
    cmd <- paste(smr_exe, "--beqtl-summary", bfile, "--query 1 --out", outfile)
    tryCatch({
      system(cmd, intern = TRUE)
    }, error = function(e) {
      rlang::abort(paste("SMR conversion failed for:", bfile))
    })
  }

  # B. Text Processing & Splitting
  process_text <- function(txt_file, out_dir, tissue, impute_maf) {
    if (!file.exists(txt_file)) rlang::abort(paste("Generated text file missing:", txt_file))

    msg <- paste0("Processing ", tissue, "...")
    message(msg)

    df <- data.table::fread(txt_file, sep = "\t", data.table = FALSE)

    # Filter
    df <- df[df$Probe %in% target_genes, ]
    if (nrow(df) == 0) {
      warning(paste("No matching genes found in", tissue))
      return(NULL)
    }

    # MAF Imputation
    if (impute_maf) {
      df <- merge(df, maf_data, by = "SNP", all.x = TRUE)
    } else if ("Freq" %in% colnames(df)) {
      # Fallback to internal freq if not imputing
      df$MAF <- ifelse(df$Freq > 0.5, 1 - df$Freq, df$Freq)
      # Check for high missingness logic can go here if strictly needed
    } else {
      rlang::abort(paste("No MAF source for", tissue))
    }

    # Rename & Select Columns
    # Map: Old -> New
    col_map <- c(SNP="SNP", Chr="CHR", BP="BP", A1="A1", A2="A2", MAF="MAF", Probe="Probe", b="BETA", SE="SE", p="P")

    # Rename matching columns
    existing <- colnames(df)
    for (old in names(col_map)) {
      if (old %in% existing) colnames(df)[colnames(df) == old] <- col_map[old]
    }

    required <- c("SNP", "CHR", "BP", "A1", "A2", "MAF", "Probe", "BETA", "SE", "P")
    if (!all(required %in% colnames(df))) rlang::abort(paste("Missing columns in", tissue))

    df <- df[, required]

    # Split & Write
    sep_dir <- file.path(out_dir, "separate_file")
    if (!dir.exists(sep_dir)) dir.create(sep_dir, recursive = TRUE)

    prefix <- if (tissue == "eQTLGEN") "eQTLGEN_" else paste0(tissue, "_")

    # Split by Probe
    lapply(split(df, df$Probe), function(sub_df) {
      data.table::fwrite(sub_df, file.path(sep_dir, paste0(prefix, sub_df$Probe[1], ".txt")), sep = "\t")
    })

    message(paste("  -> Saved gene files for", tissue))
  }

  # --- 4. Main Loop ---

  for (tissue in tissue_list) {

    # CASE 1: eQTLGen (Special Tar.gz handling)
    if (tissue == "eQTLGEN") {
      work_dir <- file.path(save_folder, "eQTLGEN")
      if (!dir.exists(work_dir)) dir.create(work_dir)

      # 1. Unpack
      tar_path <- file.path(eqtl_data_folder, "cis-eQTL-SMR_20191212.tar.gz")
      if (file.exists(tar_path)) {
        message(">>> Unpacking eQTLGen...")
        untar(tar_path, exdir = eqtl_data_folder) # Extracts to folder

        # Locate and gunzip the specific .esi/.epi/.besd files to work_dir
        # Assuming regex match for the weird long filenames
        gz_files <- list.files(eqtl_data_folder, pattern = "cis-eQTLs-full.*\\.gz$", full.names = TRUE)
        for (f in gz_files) {
          ext <- tools::file_ext(sub("\\.gz$", "", f)) # get esi/epi/besd
          R.utils::gunzip(f, destname = file.path(work_dir, paste0("eQTLGEN.", ext)), overwrite = TRUE)
        }
      }

      # 2. SMR Convert
      bfile_base <- file.path(work_dir, "eQTLGEN")
      if (!file.exists(paste0(bfile_base, ".besd"))) rlang::abort("eQTLGen binary files extraction failed.")

      run_smr(bfile_base, bfile_base) # output is eQTLGEN.txt

      # 3. Process
      process_text(paste0(bfile_base, ".txt"), work_dir, "eQTLGEN", impute_maf = FALSE) # eQTLGen has MAF/Freq usually

      # 4. Cleanup
      if (smr_file_delete) file.remove(paste0(bfile_base, c(".esi", ".epi", ".besd")))
      if (complete_file_delete) file.remove(paste0(bfile_base, ".txt"))

    } else {
      # CASE 2: GTEx
      tissue_dir <- file.path(eqtl_data_folder, tissue)
      output_dir <- file.path(save_folder, tissue)
      if (!dir.exists(output_dir)) dir.create(output_dir)

      # 1. Unzip
      zip_path <- file.path(eqtl_data_folder, paste0(tissue, ".zip"))
      message(paste(">>> Unzipping", tissue, "..."))
      unzip(zip_path, exdir = eqtl_data_folder) # Extracts to folder named 'tissue'

      # 2. SMR Convert
      bfile_base <- file.path(tissue_dir, tissue) # GTEx format: Folder/Tissue.besd
      out_base   <- file.path(tissue_dir, tissue)

      if (!file.exists(paste0(bfile_base, ".besd"))) rlang::abort(paste("Binary files missing for", tissue))

      run_smr(bfile_base, out_base)

      # 3. Process
      process_text(paste0(out_base, ".txt"), output_dir, tissue, impute_maf = maf_impute)

      # 4. Cleanup
      if (smr_file_delete) file.remove(paste0(bfile_base, c(".esi", ".epi", ".besd")))
      if (complete_file_delete) file.remove(paste0(out_base, ".txt"))
    }
  }

  message("All processing completed successfully.")
}
