#' Prepare the eQTL data used in the colocalization-based linking
#'
#' Yanglab provides a tidy dataset for GTEx v8 cis eQTLs, but this dataset lacks MAF information. Therefore, this function supplies MAF values (sourced from the “eQTL_catalogue”) to complete the GTEx summary data. In addition, this program splits the eQTL data by individual genes.
#'
#' @param eqtl_data_folder
#' @param eqtl_data_tissue
#' @param smr_format
#' @param eqtl_data_suffix
#' @param save_folder
#' @param smr_exe
#' @param extract_gene_list
#' @param smr_file_delete
#' @param complete_file_delete
#' @param maf_impute
#' @param maf_impute_file
#'
#' @return
#' @export
#'
#' @examples
eqtl_data_preparation <- function(eqtl_data_folder, eqtl_data_tissue = NULL, smr_format = TRUE, eqtl_data_suffix = NULL, save_folder = NULL,
                                  smr_exe, extract_gene_list, smr_file_delete = TRUE, complete_file_delete = TRUE, maf_impute = TRUE, maf_impute_file = NULL)
{
  suppressMessages(suppressWarnings(library(dplyr)))

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
    "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Breast_Mammary_Tissue", "Liver", "Whole_Blood", "cis-eQTL-SMR_20191212", "eqtlGEN"
  )

  if (missing(eqtl_data_folder) || !is.character(eqtl_data_folder) || length(eqtl_data_folder) != 1 || !dir.exists(eqtl_data_folder)) {
    rlang::abort("Error: 'eqtl_data_folder' must be provided and must be a valid directory path.")
  }
  eqtl_data_folder <- normalizePath(eqtl_data_folder)

  check_and_standardize_pheno_name <- function(eqtl_data_tissue) {
    if (missing(eqtl_data_tissue) || !is.character(eqtl_data_tissue) || length(eqtl_data_tissue) != 1 || nchar(eqtl_data_tissue) == 0) {
      rlang::abort("'eqtl_data_tissue' must be provided and must be a non-empty character string.")
    }
    eqtl_data_tissue <- unlist(strsplit(eqtl_data_tissue, ","))
    eqtl_data_tissue <- trimws(eqtl_data_tissue)
    if (any(nchar(eqtl_data_tissue) == 0)) {
      rlang::abort("'eqtl_data_tissue' contains empty values after splitting by commas.")
    }
    return(eqtl_data_tissue)
  }
  if(!is.null(eqtl_data_tissue)) {
    eqtl_data_tissue <- check_and_standardize_pheno_name(eqtl_data_tissue)
  }

  if (smr_format) {
    if (is.null(smr_exe)) {
      rlang::abort("Error: 'smr_exe' must be provided when 'smr_format' is TRUE.")
    }
    if (!file.exists(smr_exe)) {
      rlang::abort(paste("Error: The provided 'smr_exe' path is not valid:", smr_exe))
    }

    # check eqtl_data_suffix
    if (!is.null(eqtl_data_suffix)) {
      rlang::abort("Error: 'eqtl_data_suffix' is not allowed when 'smr_format' is TRUE. This parameter is for eQTL files already converted by SMR.")
    }

  } else {
    if (is.null(eqtl_data_suffix) || !(is.character(eqtl_data_suffix) && length(eqtl_data_suffix) == 1)) {
      rlang::abort("Error: 'eqtl_data_suffix' must be provided and must be a single string when 'smr_format' is FALSE.")
    }
  }

  zip_files <- list.files(eqtl_data_folder,
                          pattern = "\\.zip$|cis-eQTL-SMR_20191212\\.tar\\.gz$",
                          full.names = FALSE)
  zip_file <- length(zip_files) > 0
  if (zip_file) {
    if (!smr_format) {
      rlang::abort("Error: If .zip files are found, 'smr_format' must be TRUE.")
    }
    if (is.null(eqtl_data_tissue)) {
      message(paste("There are total", length(zip_files), "compressed files found in", eqtl_data_folder))
      zip_tissues <- sub("\\.zip$", "", zip_files)  # delete ".zip" suffix
      zip_tissues <- sub("\\.tar\\.gz$", "", zip_tissues)  # delete ".tar.gz" suffix
      tissue_list <- zip_tissues
    } else {
      zip_tissues <- sub("\\.zip$", "", zip_files)  # delete ".zip" suffix
      zip_tissues <- sub("\\.tar\\.gz$", "", zip_tissues)  # delete ".tar.gz" suffix
      invalid_tissues <- eqtl_data_tissue[!eqtl_data_tissue %in% zip_tissues]
      if (length(invalid_tissues) > 0) {
        rlang::abort(paste("Error: The following tissues are not found as compressed format in the folder:",
                   paste(invalid_tissues, collapse = ", ")))
      }
      tissue_list <- eqtl_data_tissue
    }
    invalid_tissues <- tissue_list[!tissue_list %in% valid_tissue_names]
    if (length(invalid_tissues) > 0) {
      warning(paste0("warning: The following tissues are not in the GTeX-V1/eqtlGEN-V8 valid tissue list:", " ",
                     paste(invalid_tissues, collapse = ", "),". please check."))
    }
  } else {
    if (is.null(eqtl_data_tissue)) {
      tissue_list <- list.dirs(eqtl_data_folder, recursive = FALSE, full.names = FALSE)
      invalid_tissues <- tissue_list[!tissue_list %in% valid_tissue_names]
      if (length(invalid_tissues) > 0) {
        warning(paste0("warning: The following tissues are not in the GTeX-V1/eqtlGEN-V8 valid tissue list:", " ",
                            paste(invalid_tissues, collapse = ", "),". please check."))
      }
      if (smr_format) {
        valid_tissues <- tissue_list[sapply(tissue_list, function(x) {
          if (x != "cis-eQTL-SMR_20191212") {
            all(file.exists(file.path(eqtl_data_folder, x, paste0(x, c(".esi", ".epi", ".besd")))))
          } else {
            esi_exists <- file.exists(file.path(eqtl_data_folder, x,
                                                "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi")) ||
              file.exists(file.path(eqtl_data_folder, x,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi.gz"))

            epi_exists <- file.exists(file.path(eqtl_data_folder, x,
                                                "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi")) ||
              file.exists(file.path(eqtl_data_folder, x,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi.gz"))

            besd_exists <- file.exists(file.path(eqtl_data_folder, x,
                                                 "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd")) ||
              file.exists(file.path(eqtl_data_folder, x,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd.gz"))

            esi_exists && epi_exists && besd_exists
          }
        })]

        message(paste("There are total", length(valid_tissues), "valid tissues in smr format in", eqtl_data_folder))
        tissue_list <- valid_tissues
      } else {
        valid_tissues <- tissue_list[sapply(tissue_list, function(x) {
          if (x != "cis-eQTL-SMR_20191212") {
            file.exists(file.path(eqtl_data_folder, x, paste0(x, eqtl_data_suffix)))
          } else {
            esi_exists <- file.exists(file.path(eqtl_data_folder, x,
                                                "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi")) ||
              file.exists(file.path(eqtl_data_folder, x,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi.gz"))

            epi_exists <- file.exists(file.path(eqtl_data_folder, x,
                                                "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi")) ||
              file.exists(file.path(eqtl_data_folder, x,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi.gz"))

            besd_exists <- file.exists(file.path(eqtl_data_folder, x,
                                                 "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd")) ||
              file.exists(file.path(eqtl_data_folder, x,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd.gz"))

            esi_exists && epi_exists && besd_exists
          }
        })]
        tissue_list <- valid_tissues
        message(paste("There are total", length(valid_tissues), "tissues in txt format in", eqtl_data_folder))
      }
    } else {
      invalid_tissues <- eqtl_data_tissue[!eqtl_data_tissue %in% valid_tissue_names]
      if (length(invalid_tissues) > 0) {
        warning(paste0("warning: The following tissues are not in the GTeX-V1/eqtlGEN-V8 valid tissue list:", " ",
                            paste(invalid_tissues, collapse = ", "),". please check."))
      }
      if (smr_format) {
        for (tissue in eqtl_data_tissue) {
          if (tissue != "cis-eQTL-SMR_20191212") {
            if (!(file.exists(file.path(eqtl_data_folder, tissue, paste0(tissue, ".esi"))) &&
                  file.exists(file.path(eqtl_data_folder, tissue, paste0(tissue, ".epi"))) &&
                  file.exists(file.path(eqtl_data_folder, tissue, paste0(tissue, ".besd"))))) {
              rlang::abort(paste("Error: The tissue", tissue, "is missing one or more of the required SMR format files (.esi, .epi, .besd)."))
            }
          } else {
            esi_exists <- file.exists(file.path(eqtl_data_folder, tissue,
                                                "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi")) ||
              file.exists(file.path(eqtl_data_folder, tissue,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi.gz"))

            epi_exists <- file.exists(file.path(eqtl_data_folder, tissue,
                                                "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi")) ||
              file.exists(file.path(eqtl_data_folder, tissue,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi.gz"))

            besd_exists <- file.exists(file.path(eqtl_data_folder, tissue,
                                                 "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd")) ||
              file.exists(file.path(eqtl_data_folder, tissue,
                                    "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd.gz"))

            if (!(esi_exists && epi_exists && besd_exists)) {
              rlang::abort(paste("Error: The tissue", tissue, "is missing one or more of the required files for SMR format (either .esi, .epi, .besd or their .gz counterparts)."))
            }
          }
        }
      } else {
        for (tissue in eqtl_data_tissue) {
          if (!file.exists(file.path(eqtl_data_folder, tissue, paste0(tissue, eqtl_data_suffix)))) {
            rlang::abort(paste("Error: The tissue", tissue, "is missing the required file with suffix", eqtl_data_suffix))
          }
        }
      }
      tissue_list <- eqtl_data_tissue
    }
  }

  ##check save_folder
  if (!is.null(save_folder)) {
    if (missing(save_folder) || !is.character(save_folder) || length(save_folder) != 1 || !dir.exists(save_folder)) {
      rlang::abort("'save_folder' must be NULL(same as eqtl_data_folder) or a valid directory path.")
    }
    if (!file.access(save_folder, 2) == 0) {
      rlang::abort("'save_folder' must be a directory where the user has write permissions.")
    }
    save_folder <- normalizePath(save_folder)
  } else {
    save_folder <- eqtl_data_folder
  }


  if (missing(extract_gene_list) || !is.character(extract_gene_list) || length(extract_gene_list) != 1 || !file.exists(extract_gene_list)) {
    rlang::abort("Error: 'eqtl_data_folder' must be provided and must be a valid file")
  }
  extract_gene_list <- data.table::fread(extract_gene_list, header = FALSE, data.table = FALSE)
  if (!(is.data.frame(extract_gene_list) || is.matrix(extract_gene_list)) || ncol(extract_gene_list) != 1) {
    rlang::abort("'extract_gene_list' must be a data frame or matrix with single column.")
  }


  if (!is.logical(smr_file_delete) || length(smr_file_delete) != 1) {
    rlang::abort("'smr_file_delete' must be either TRUE or FALSE.")
  }

  if (!is.logical(complete_file_delete) || length(complete_file_delete) != 1) {
    rlang::abort("'complete_file_delete' must be either TRUE or FALSE.")
  }

  if (!is.logical(maf_impute) || length(maf_impute) != 1) {
    rlang::abort("'maf_impute' must be either TRUE or FALSE.")
  }

  if (!maf_impute && !is.null(maf_impute_file)) {
    rlang::abort("Error: When maf_impute is FALSE, maf_impute_file must be NULL.")
  }
  if (maf_impute && is.null(maf_impute_file)) {
    maf_data <- data.table::fread(system.file("extdata", "maf_impute.csv.gz", package = "mds2g"), data.table = FALSE, stringsAsFactors = FALSE, header=T)
  } else {
    if (!is.character(maf_impute_file) || length(maf_impute_file) != 1 || !file.exists(maf_impute_file)) {
      rlang::abort("Error: maf_impute_file must be a valid character string pointing to an existing file.")
    }
    maf_data <- data.table::fread(maf_impute_file, data.table = FALSE, stringsAsFactors = FALSE, header=T)
    if (!(is.data.frame(maf_data) || is.matrix(maf_data)) || ncol(maf_data) < 2) {
      rlang::abort("Error: The specified maf_impute_file must be a data.frame or matrix with at least two columns(refer to rsid and maf).")
    }
    if (!all(c("SNP", "MAF") %in% colnames(maf_data))) {
      rlang::abort("Error: The specified maf_impute_file must contain 'rsid' and 'maf' columns.")
    }
    maf_data <- maf_data[, c("SNP", "MAF")]
  }



  data_preparation <- function(eqtl_data_folder, tissue, smr_format, eqtl_data_suffix, save_folder,
                               smr_exe, zip_file, extract_gene_list, smr_file_delete, complete_file_delete, maf_impute, maf_impute_file)
  {
    if (zip_file) {
      zip_path <- file.path(eqtl_data_folder, paste0(tissue, ".zip"))
      if (file.exists(zip_path)) {
        unzip(zip_path, exdir = eqtl_data_folder)
        file.remove(zip_path)
        message(paste("Unzipped", zip_path, "to", eqtl_data_folder))
      } else {
        rlang::abort(paste("Error: Zip file", zip_path, "does not exist."))
      }
    }

    if (smr_format) {
      smr_cmd <- paste(smr_exe, "--beqtl-summary", file.path(eqtl_data_folder, tissue, tissue),
                       "--query 1 --out", file.path(eqtl_data_folder, tissue, tissue))
      result <- tryCatch({
        system(smr_cmd, intern = TRUE)
      }, error = function(e) {
        txt_file <- file.path(eqtl_data_folder, tissue, paste0(tissue, ".txt"))
        if (file.exists(txt_file)) file.remove(txt_file)
        rlang::abort(paste("Error: SMR command failed for", tissue))
      })
      message(paste("converted SMR-formatted file to", file.path(eqtl_data_folder, tissue, paste0(tissue, ".txt"))))

      if (smr_file_delete) {
        file.remove(file.path(eqtl_data_folder, tissue, paste0(tissue, ".esi")),
                    file.path(eqtl_data_folder, tissue, paste0(tissue, ".epi")),
                    file.path(eqtl_data_folder, tissue, paste0(tissue, ".besd")))
      }
    }

    if (smr_format) {
      data_file <- file.path(eqtl_data_folder, tissue, paste0(tissue, ".txt"))
    } else {
      data_file <- file.path(eqtl_data_folder, tissue, paste0(tissue, eqtl_data_suffix))
    }

    if (!file.exists(data_file)) rlang::abort(paste("Error: Data file", data_file, "does not exist."))

    eqtl_data <- data.table::fread(data_file, sep = "\t", data.table = FALSE, stringsAsFactors = FALSE)

    eqtl_data <- eqtl_data[eqtl_data$Probe %in% extract_gene_list$V1, ]

    if (maf_impute) {
      eqtl_data <- merge(eqtl_data, maf_data, by = "SNP", all.x = TRUE)
    } else {
      freq_na_ratio <- sum(is.na(eqtl_data$Freq)) / nrow(eqtl_data)
      if (freq_na_ratio <= 0.8) {
        eqtl_data$MAF <- ifelse(eqtl_data$Freq > 0.5, 1 - eqtl_data$Freq, eqtl_data$Freq)
      } else if (freq_na_ratio > 0.8) {
        rlang::abort(paste0("Error: More than 80% of the Freq column in ", tissue, " is NA, please consider setting 'maf_impute' to TRUE."))
      }
    }

    eqtl_data <- eqtl_data %>%
      dplyr::select(-Probe_Chr, -Probe_bp, -Gene, -Orientation, -Freq) %>%
      dplyr::select(SNP, Chr, BP, A1, A2, MAF, Probe, b, SE, p) %>%
      dplyr::rename(
        SNP = SNP,
        CHR = Chr,
        BP = BP,
        A1 = A1,
        A2 = A2,
        MAF = MAF,
        Probe = Probe,
        BETA = b,
        SE = SE,
        P = p
      )

    separate_folder <- file.path(save_folder, tissue, "separate_file")
    if (!dir.exists(separate_folder)) dir.create(separate_folder, recursive = TRUE)

    unique_probes <- unique(eqtl_data$Probe)

    for (probe in unique_probes) {
      probe_data <- eqtl_data[eqtl_data$Probe == probe, ]
      probe_file <- file.path(separate_folder, paste0(tissue, "_", probe, ".txt"))
      data.table::fwrite(probe_data, probe_file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }

    if (complete_file_delete) {
      file.remove(data_file)
    }

    message(paste("Data preparation complete for", tissue))
  }




  data_preparation_eqtlGEN <- function(eqtl_data_folder, tissue, smr_format, eqtl_data_suffix, save_folder,
                               smr_exe, zip_file, extract_gene_list, smr_file_delete, complete_file_delete, maf_impute, maf_impute_file)
  {
    eqtlGEN <- "eqtlGEN"
    if (zip_file) {
      tar_path <- file.path(eqtl_data_folder, "cis-eQTL-SMR_20191212.tar.gz")
      if (file.exists(tar_path)) {
        untar(tar_path, exdir = eqtl_data_folder)
        file.remove(tar_path)
        message(paste("Unzipped and removed", tar_path))
        esi_gz_path <- file.path(eqtl_data_folder,
                                 "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi.gz")
        epi_gz_path <- file.path(eqtl_data_folder,
                                 "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi.gz")
        besd_gz_path <- file.path(eqtl_data_folder,
                                  "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd.gz")
        eqtlGEN_folder <- file.path(save_folder, eqtlGEN)
        if (!dir.exists(eqtlGEN_folder)) dir.create(eqtlGEN_folder, recursive = TRUE)
        if (file.exists(esi_gz_path)) {
          R.utils::gunzip(esi_gz_path, destname = file.path(eqtlGEN_folder, "eqtlGEN.esi"), overwrite = TRUE)
        }
        if (file.exists(epi_gz_path)) {
          R.utils::gunzip(epi_gz_path, destname = file.path(eqtlGEN_folder, "eqtlGEN.epi"), overwrite = TRUE)
        }
        if (file.exists(besd_gz_path)) {
          R.utils::gunzip(besd_gz_path, destname = file.path(eqtlGEN_folder, "eqtlGEN.besd"), overwrite = TRUE)
        }

        message("Successfully unzipped and renamed files for cis-eQTL-SMR_20191212.")
      } else {
        rlang::abort(paste("Error: Tar file", tar_path, "does not exist."))
      }
    }

    if (!zip_file && smr_format) {
      smr_folder <- file.path(eqtl_data_folder, tissue)
      smr_files_gz <- list.files(smr_folder, pattern = "\\.gz$", full.names = TRUE)
      eqtlGEN_folder <- file.path(save_folder, eqtlGEN)
      if (!dir.exists(eqtlGEN_folder)) dir.create(eqtlGEN_folder, recursive = TRUE)
      if (length(smr_files_gz) > 0) {
        esi_gz <- grep("\\.esi\\.gz$", smr_files_gz, value = TRUE)
        epi_gz <- grep("\\.epi\\.gz$", smr_files_gz, value = TRUE)
        besd_gz <- grep("\\.besd\\.gz$", smr_files_gz, value = TRUE)

        if (length(esi_gz) == 0 || length(epi_gz) == 0 || length(besd_gz) == 0) {
          rlang::abort("Error: Missing one or more of the required .gz files (.esi.gz, .epi.gz, .besd.gz)")
        }

        R.utils::gunzip(esi_gz, destname = file.path(eqtlGEN_folder, "eqtlGEN.esi"), overwrite = TRUE)
        R.utils::gunzip(epi_gz, destname = file.path(eqtlGEN_folder, "eqtlGEN.epi"), overwrite = TRUE)
        R.utils::gunzip(besd_gz, destname = file.path(eqtlGEN_folder, "eqtlGEN.besd"), overwrite = TRUE)

      } else {
        esi_file <- file.path(smr_folder, "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi")
        epi_file <- file.path(smr_folder, "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi")
        besd_file <- file.path(smr_folder, "cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd")

        if (!(file.exists(esi_file) && file.exists(epi_file) && file.exists(besd_file))) {
          rlang::abort("Error: Missing one or more of the required SMR files (.esi, .epi, .besd)")
        }

        file.rename(esi_file, file.path(eqtlGEN_folder, "eqtlGEN.esi"))
        file.rename(epi_file, file.path(eqtlGEN_folder, "eqtlGEN.epi"))
        file.rename(besd_file, file.path(eqtlGEN_folder, "eqtlGEN.besd"))
      }
    }

    smr_cmd <- paste(smr_exe, "--beqtl-summary", file.path(eqtlGEN_folder, "eqtlGEN"),
                     "--query 1 --out", file.path(eqtlGEN_folder, "eqtlGEN"))
    result <- tryCatch({
      system(smr_cmd, intern = TRUE)
    }, error = function(e) {
      rlang::abort(paste("Error: SMR command failed for", tissue))
    })
    if (smr_file_delete) {
      file.remove(file.path(eqtlGEN_folder, "eqtlGEN.esi"),
                  file.path(eqtlGEN_folder, "eqtlGEN.epi"),
                  file.path(eqtlGEN_folder, "eqtlGEN.besd"))
    }
    message(paste("converted SMR-formatted file to", file.path(eqtlGEN_folder, paste0(eqtlGEN, ".txt"))))

    if (!zip_file && !smr_format) {
      data_file <- file.path(eqtl_data_folder, tissue, paste0(tissue, eqtl_data_suffix))

      if (!file.exists(data_file)) {
        rlang::abort(paste("Error: Data file", data_file, "does not exist."))
      }

      eqtlGEN_folder <- file.path(eqtl_data_folder, eqtlGEN)
      if (!dir.exists(eqtlGEN_folder)) dir.create(eqtlGEN_folder, recursive = TRUE)

      file.rename(data_file, file.path(eqtlGEN_folder, "eqtlGEN.txt"))

      message("Successfully moved and renamed data file for", tissue)
    }

    eqtlGEN_file <- file.path(eqtl_data_folder, eqtlGEN, "eqtlGEN.txt")
    if (!file.exists(eqtlGEN_file)) rlang::abort(paste("Error: Data file", eqtlGEN_file, "does not exist."))

    eqtl_data <- data.table::fread(eqtlGEN_file, sep = "\t", data.table = FALSE, stringsAsFactors = FALSE)


    eqtl_data <- eqtl_data[eqtl_data$Probe %in% extract_gene_list$V1, ]

    if (maf_impute) {
      eqtl_data <- merge(eqtl_data, maf_data, by = "SNP", all.x = TRUE)
    } else {
      freq_na_ratio <- sum(is.na(eqtl_data$Freq)) / nrow(eqtl_data)
      if (freq_na_ratio <= 0.8) {
        eqtl_data$MAF <- ifelse(eqtl_data$Freq > 0.5, 1 - eqtl_data$Freq, eqtl_data$Freq)
      } else if (freq_na_ratio > 0.8) {
        rlang::abort(paste0("Error: More than 80% of the Freq column in ", tissue, " is NA, please consider setting 'maf_impute' to TRUE."))
      }
    }

    eqtl_data <- eqtl_data %>%
      dplyr::select(-Probe_Chr, -Probe_bp, -Gene, -Orientation, -Freq) %>%
      dplyr::select(SNP, Chr, BP, A1, A2, MAF, Probe, b, SE, p) %>%
      dplyr::rename(
        SNP = SNP,
        CHR = Chr,
        BP = BP,
        A1 = A1,
        A2 = A2,
        MAF = MAF,
        Probe = Probe,
        BETA = b,
        SE = SE,
        P = p
      )

    separate_folder <- file.path(save_folder, eqtlGEN, "separate_file")
    if (!dir.exists(separate_folder)) dir.create(separate_folder, recursive = TRUE)

    unique_probes <- unique(eqtl_data$Probe)

    for (probe in unique_probes) {
      probe_data <- eqtl_data[eqtl_data$Probe == probe, ]
      probe_file <- file.path(separate_folder, paste0("eqtlGEN_", probe, ".txt"))
      data.table::fwrite(probe_data, probe_file, sep = "\t", row.names = FALSE, col.names = TRUE)
    }

    if (complete_file_delete) {
      file.remove(eqtlGEN_file)
    }

    message(paste("Data preparation complete for", tissue))
  }

  for (tissue in tissue_list) {

    if (tissue == "cis-eQTL-SMR_20191212") {
      data_preparation_eqtlGEN(eqtl_data_folder = eqtl_data_folder,
                               tissue = tissue,
                               smr_format = smr_format,
                               eqtl_data_suffix = eqtl_data_suffix,
                               smr_exe = smr_exe,
                               save_folder = save_folder,
                               zip_file = zip_file,
                               extract_gene_list = extract_gene_list,
                               smr_file_delete = smr_file_delete,
                               complete_file_delete = complete_file_delete,
                               maf_impute = maf_impute,
                               maf_impute_file = maf_impute_file)

    } else {
      data_preparation(eqtl_data_folder = eqtl_data_folder,
                       tissue = tissue,
                       smr_format = smr_format,
                       eqtl_data_suffix = eqtl_data_suffix,
                       smr_exe = smr_exe,
                       save_folder = save_folder,
                       zip_file = zip_file,
                       extract_gene_list = extract_gene_list,
                       smr_file_delete = smr_file_delete,
                       complete_file_delete = complete_file_delete,
                       maf_impute = maf_impute,
                       maf_impute_file = maf_impute_file)
    }
  }
}



