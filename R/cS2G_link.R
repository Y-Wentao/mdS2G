#' Combined S2G Strategy
#'
#' Use locus-based gene linking by combined SNP-to-GENE (cS2G) strategy
#'
#' @param pheno_name The pheno or trait include in your study, for example "height"; If you have two or more pheno, please separate them with commas, for example "height,BMI".
#' @param SNP_data The path point to your causal SNP or lead SNP (it is not recommended yet) result, but you should never use them together. The data should have each line representing a SNP, and it should contain at least CHR, SNP, BP, PHENO (if multiple phenotypes are provided in pheno_name), and PIPS (if they are causal SNPs).
#' @param s2g_folder The path point to 'S2G_bed' mapping file. This file should be compressed and contain "Exon", "Promoter", "finemappedciseQTLs_eQTLGen", finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood" seven subfolders. The specific files within these subfolders should be uncompressed
#' @param save_folder The path where you want to store the result file.
#' @param keep_intermediate_file Whether you want to store intermediate files, default to False.
#' @param intermediate_file_folder If you want to store intermediate files, you should provide a specific storage path. It's default to NULL
#' @param gene_remove_list We suggest deleting non coding genes to clarify the results, as we have included a list of genes that should be removed. If you do not want to exclude, please select False; If you wish to add an additional list, please provide the specific path of the single column text without a header. It's default to NULL.
#' @param other_weights We used the optimal weights provided in the Gazal manuscript as the integrated weights for the seven strategies. If you have any other weights you would like to use, please provide the specific file. This data file should not have a header. The first column corresponds to the names of the seven strategies used (Please refer to the description of the s2d_folder parameter above for details), the second column should correspond to specific weights. That is to say, this file should have seven rows and two columns. It's default to NULL.
#' @param other_precisions We used the precision of these seven strategies from the Gazal manuscript in the validation set as part of the confidence score calculation. If you have any other data you would like to use, please enter it. The specific requirements of the document are the same as' other weights'. It's default to NULL.
#' @param causal_SNP Whether the SNPs included in 'SNP_data' are causal SNPs, default to TRUE.
#' @param assumed_PIP If your SNP does not come from fine mapping, you need to provide a hypothetical posterior probability to calculate the confidence score. We recommend using a conservative data (default is 0.5), although this may require support from the other two types of evidence for all linking results, it can effectively reduce false positive results.
#' @param gene_annot_file In order to be compatible with the other two types of strategies, we need to convert the gene names obtained from s2g into ENSGID. We have provided the corresponding data, but if you want to use another list, please specify it. It's default to NULL.
#' @param keep_unmatched_gene Whether to retain genes that failed to convert ENSGID, default to False
#'
#' @return A data.frame named pheno_cS2G.csv
#' @export

cs2g_link <- function(pheno_name, SNP_data, s2g_folder, save_folder, keep_intermediate_file = FALSE, intermediate_file_folder = NULL,
         gene_remove_list = NULL, other_weights = FALSE, other_precisions = FALSE, causal_SNP = TRUE, assumed_PIP = NULL,
         gene_annot_file = NULL, keep_unmatched_gene = FALSE)
{
  suppressMessages(suppressWarnings(library(dplyr)))
  # check pheno_name

    check_and_standardize_pheno_name <- function(pheno_name) {
    if (missing(pheno_name) || !is.character(pheno_name) || length(pheno_name) != 1 || nchar(pheno_name) == 0) {
      rlang::abort("'pheno_name' must be provided and must be a non-empty character string.")
    }
    pheno_name <- unlist(strsplit(pheno_name, ","))
    pheno_name <- trimws(pheno_name)
    if (any(nchar(pheno_name) == 0)) {
      rlang::abort("'pheno_name' contains empty values after splitting by commas.")
    }
    return(pheno_name)
  }

  # Check causal_SNP
  if (!is.logical(causal_SNP) || length(causal_SNP) != 1) {
    rlang::abort("'causal_SNP' must be TRUE or FALSE.")
  }

  # check SNP_data
  check_and_standardize_SNP_data <- function(SNP_data, causal_SNP, pheno_name) {
    if (missing(SNP_data)) {
      rlang::abort("'SNP_data' must be provided as a file path.")
    }
    SNP_data <- read.csv(SNP_data, header=T)
    if (!is.data.frame(SNP_data) && !is.matrix(SNP_data)) {
      rlang::abort("'SNP_data' must be a valid file path or a data frame/matrix.")
    }
    colnames_lower <- tolower(colnames(SNP_data))
    required_columns <- list(
      CHR = c("chr", "chrom"),
      SNP = "snp",
      BP = c("bp", "pos")
    )
    chr_col <- colnames(SNP_data)[colnames_lower %in% tolower(required_columns$CHR)]
    snp_col <- colnames(SNP_data)[colnames_lower %in% tolower(required_columns$SNP)]
    bp_col <- colnames(SNP_data)[colnames_lower %in% tolower(required_columns$BP)]
    if (length(chr_col) == 0 || length(snp_col) == 0 || length(bp_col) == 0) {
      rlang::abort("SNP_data must contain columns for 'CHR', 'SNP', and 'BP' (or equivalents).")
    }

    pip_col <- NULL
    if (causal_SNP) {
      pip_candidates <- c("pip", "pips", "pp")
      pip_col <- colnames(SNP_data)[colnames_lower %in% pip_candidates]
      if (length(pip_col) == 0) {
        rlang::abort("For causal SNP analysis, SNP_data must contain a 'PIP' (or equivalent) column.")
      }
    }

    # only length of pheno_name > 1, check PHENO column
    pheno_col <- NULL
    if (length(pheno_name) > 1) {
      pheno_candidates <- c("pheno", "phenotype", "trait")
      pheno_col <- colnames(SNP_data)[colnames_lower %in% pheno_candidates]
      if (length(pheno_col) == 0) {
        rlang::abort("SNP_data must contain a 'PHENO' (or equivalent) column when multiple 'pheno_name' are provided.")
      }

      # standardize the column name
      selected_cols <- c(chr_col, snp_col, bp_col, if (!is.null(pip_col)) pip_col, pheno_col)
      SNP_data <- SNP_data[, selected_cols]
      colnames(SNP_data) <- c("CHR", "SNP", "BP", if (!is.null(pip_col)) "PIP", "PHENO")

      # check whether the pheno_name is consistent with PHENO column of SNP_data
      unique_pheno_in_data <- unique(SNP_data$PHENO)
      if (!all(unique(pheno_name) %in% unique_pheno_in_data)) {
        unmatched_pheno <- setdiff(pheno_name, unique_pheno_in_data)
        rlang::abort(paste("The following 'pheno_name' entries are not found in 'SNP_data':", paste(unmatched_pheno, collapse = ", ")))
      }
    } else {
      # PHENO column is not need when pheno_name is single
      selected_cols <- c(chr_col, snp_col, bp_col, if (!is.null(pip_col)) pip_col)
      SNP_data <- SNP_data[, selected_cols]
      colnames(SNP_data) <- c("CHR", "SNP", "BP", if (!is.null(pip_col)) "PIP")
    }

    # standardize the CHR column
    SNP_data$CHR <- gsub("chr", "", SNP_data$CHR, ignore.case = TRUE)  # remove 'chr' prefix
    SNP_data$CHR <- as.numeric(trimws(SNP_data$CHR))  # convert 'CHR' column to numeric
    SNP_data <- SNP_data[SNP_data$CHR %in% 1:22, ]  # keep autosomal chromosome

    return(SNP_data)
  }

  pheno_name <- check_and_standardize_pheno_name(pheno_name)
  SNP_data <- check_and_standardize_SNP_data(SNP_data, causal_SNP = causal_SNP, pheno_name = pheno_name)

  # check cs2g_folder
  if (missing(s2g_folder) || !is.character(s2g_folder) || length(s2g_folder) != 1 || !dir.exists(s2g_folder)) {
    rlang::abort("'s2g_folder' must be provided and must be a valid directory path.")
  }
  s2g_folder <- normalizePath(s2g_folder)
  required_subfolders <- c("Exon", "Promoter", "finemappedciseQTLs_eQTLGen",
                           "finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood")
  missing_subfolders <- required_subfolders[!dir.exists(file.path(s2g_folder, required_subfolders))]
  if (length(missing_subfolders) > 0) {
    rlang::abort(paste("cS2G annotation files are not as expected, please check:", paste(missing_subfolders, collapse = ", ")))
  }

  # check cs2g_file
  for (subfolder in required_subfolders) {
    for (chr in 1:22) {
      bed_file <- file.path(s2g_folder, subfolder, paste0("chr", chr, ".bed.gz"))
      if (!file.exists(bed_file)) {
        rlang::abort(paste0("Missing file: ", bed_file, ". cS2G annotation files are not as expected, please check."))
      }
    }
  }

  # check save_folder
  if (missing(save_folder) || !is.character(save_folder) || length(save_folder) != 1 || !dir.exists(save_folder)) {
    rlang::abort("'save_folder' must be provided and must be a valid directory path.")
  }
  if (!file.access(save_folder, 2) == 0) {
    rlang::abort("'save_folder' must be a directory where the user has write permissions.")
  }
  save_folder <- normalizePath(save_folder)

  # check keep_intermediate_file
  if (!is.logical(keep_intermediate_file) || length(keep_intermediate_file) != 1) {
    rlang::abort("'keep_intermediate_file' must be TRUE or FALSE!")
  }
  # If keep_intermediate_file is TRUE, check intermediate_file_folder
  if (keep_intermediate_file) {
    if (is.null(intermediate_file_folder)) {
      rlang::abort("'intermediate_file_folder' must be provided when 'keep_intermediate_file' is TRUE.")
    }
    if (!file.access(intermediate_file_folder, 2) == 0) {
      rlang::abort("'intermediate_file_folder' must be a directory where the user has write permissions.")
    }
    intermediate_file_folder <- normalizePath(intermediate_file_folder)
  }
  if(!is.null(intermediate_file_folder)){
    if(!keep_intermediate_file){
      rlang::abort("'Keep_intermediate_file' must be TRUE when 'intermediate_file_folder' is provided.")
    }
  }

  # Check gene_remove_list
  if (isTRUE(gene_remove_list)) {
    rlang::abort("'gene_remove_list' must be NULL, FALSE or clear gene list.")
  }
  if (!is.null(gene_remove_list) && is.character(gene_remove_list)) {
    gene_remove_list <- data.table::fread(gene_remove_list, header = F, data.table=FALSE)
    if (!(is.data.frame(gene_remove_list) || is.matrix(gene_remove_list))) {
      rlang::abort("'gene_remove_list' must be a data frame or matrix.")
    }
    gene_remove_list <- gene_remove_list %>% dplyr::pull(1)
  }

  # Check other_weights
  if (isTRUE(other_weights)) {
    rlang::abort("'other_weights' must be NULL, FALSE or clear weights list.")
  }
  if (!identical(other_weights, FALSE)) {
    other_weights_data <- data.table::fread(other_weights, header = F, data.table=FALSE)
    if (nrow(other_weights_data) < 7 || ncol(other_weights_data) != 2) {
      rlang::abort("'other_weights' must be a data frame or matrix with 2 columns and at least 7 rows.")
    }
    weights_vector <- setNames(as.numeric(other_weights_data[[2]]), other_weights_data[[1]])
    required_weights <- c("Exon", "Promoter", "finemappedciseQTLs_eQTLGen", "finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood")
    if (!all(required_weights %in% names(weights_vector))) {
      rlang::abort("The vector names for 'other_weights' must include Exon, Promoter, finemappedciseQTLs_eQTLGen, finemappedciseQTLs_GTeX, EpiMap, ABC, Ciceroblood, please check")
    }
    if (!is.numeric(weights_vector)) {
      rlang::abort("'other_weights' must contain numeric values.")
    }
  }

  # Check other_precisions
  if (isTRUE(other_precisions)) {
    rlang::abort("'other_precisions' must be NULL, FALSE or clear precisions list.")
  }

  if (!identical(other_precisions, FALSE)) {
    other_precisions_data <- data.table::fread(other_precisions, header = F, data.table=FALSE)
    if (nrow(other_precisions_data) < 7 || ncol(other_precisions_data) < 2) {
      rlang::abort("'other_precisions' must be a data frame or matrix with at least 7 rows and 2 columns.")
    }
    precisions_vector <- setNames(as.numeric(other_precisions_data[[2]]), other_precisions_data[[1]])
    if (!all(required_weights %in% names(precisions_vector))) {
      rlang::abort("The vector names for 'other_precisions' must include Exon, Promoter, finemappedciseQTLs_eQTLGen, finemappedciseQTLs_GTeX, EpiMap, ABC, Ciceroblood, please check")
    }
    if (!is.numeric(precisions_vector)) {
      rlang::abort("'other_precisions' must contain numeric values.")
    }
  }

  # If causal_SNP is FALSE, check assumed_PIP
  if (!causal_SNP) {
    if (is.null(assumed_PIP)) {
      warning("assumed_PIP not provided, defaulting to 0.5. Please check your input.")
      assumed_PIP <- 0.5 # default value
    }
    if (!is.numeric(assumed_PIP)) {
      rlang::abort("`assumed_PIP` must be a numeric value.")
    }
  } else {
    if (!is.null(assumed_PIP)) {
      rlang::abort("`causal_SNP` must be FALSE when `assumed_PIP` is provided.")
    }
  }

  # check GENE_annot file
  if(!is.null(gene_annot_file)){
    if(!is.character(gene_annot_file) || length(gene_annot_file) != 1){
      rlang::abort("'gene_annot_file' is invalid, please check.")
    } else {
      data.table::fread(gene_annot_file, header=T, stringsAsFactors = F, data.table=F)
      if (!(is.data.frame(gene_annot_file) || is.matrix(gene_annot_file)) || ncol(other_weights_data) != 2) {
        rlang::abort("'gene_annot_file' must be a data frame or matrix with at least 2 columns(refer to ENSGID and NAME)")
      }
      required_columns <- c("ENSGID", "NAME")
      if (!all(required_columns %in% colnames(gene_annot_data))) {
        rlang::abort("'gene_annot_file' must contain the columns 'ENSGID' and 'NAME'.")
      }
      gene_annot_data <- unique(gene_annot_data[, required_columns, drop = FALSE])
    }
  } else {
    gene_annot_file <- data.table::fread(system.file("extdata", "GENE_annot.csv.gz", package = "mds2g"), header=T, stringsAsFactors = F, data.table=F) %>%
      dplyr::select(ENSGID, NAME) %>% unique()
  }

  # check keep_unmatched_gene
  if (!is.logical(keep_unmatched_gene) || length(keep_unmatched_gene) != 1) {
    rlang::abort("'keep_unmatched_gene' must be TRUE or FALSE.")
  }


  link <- function (pheno_name, SNP_data, s2g_folder, save_folder, keep_intermediate_file = FALSE, intermediate_file_folder =FALSE,
                         gene_remove_list = NULL,other_weights = FALSE, other_precisions = FALSE, causal_SNP = TRUE, assumed_PIP = NULL,
                         gene_annot_file = NULL, keep_unmatched_gene = FALSE)
  {
    ####### cs2g step1:map SNP to GENE
    pheno <- pheno_name
    dir_list <- c("Exon", "Promoter", "finemappedciseQTLs_eQTLGen", "finemappedciseQTLs_GTeX", "EpiMap", "ABC", "Ciceroblood")
    # seven sub strategy loop
    for(dir in dir_list) {
      data_combination <- data.frame(SNP=character(), GENE=character(), VALUE=numeric(), STRATEGY=character(), stringsAsFactors=FALSE)
      chr_data_cache <- list()
      for (i in 1:nrow(SNP_data)) {
        # get the current file path
        chr_gz_file <- file.path(s2g_folder, dir, paste0("chr", SNP_data[i, 1], ".bed.gz"))

        # If this file has not been read yet, read it
        if (!chr_gz_file %in% names(chr_data_cache)) {
          if (file.exists(chr_gz_file)) {
            chr_data_cache[[chr_gz_file]] <- data.table::fread(chr_gz_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE)
          } else {
            rlang::abort(paste0("file ", chr_gz_file, " is missing!"))
          }
        }
        chr_data <- chr_data_cache[[chr_gz_file]]
        pos <- SNP_data[i, 3]
        in_interval <- chr_data[pos >= chr_data[[2]] & pos <= chr_data[[3]], ]
        if (nrow(in_interval) == 0) {
          # if there are no valid intervals, add NA
          data_combination <- rbind(data_combination, c(SNP_data[i, 2], NA, NA, NA))
        } else {
          # add to the data_combination
          for (j in 1:nrow(in_interval)) {
            data_combination <- rbind(data_combination, c(SNP_data[i, 2], in_interval[j, 4], in_interval[j, 5], dir))
          }
        }
      }
      # change the column name
      colnames(data_combination) <- c("SNP", "GENE", "VALUE", "STRATEGY")
      data_combination$VALUE <- as.numeric(data_combination$VALUE)
      data_combination <- data_combination[!duplicated(data_combination),]
      replace_dict <- c("MARC1" = "MTARC1",
                        "SEPT9" = "SEPTIN9",
                        "SEPT7" = "SEPTIN7",
                        "NAGR2" = "ICE2",
                        "LINC00521" = "CCDC197",
                        "FAM211B" = "LRRC75B",
                        "PCDP1" = "CFAP221",
                        "KIAA1432" = "RIC1",
                        "LINC00493" = "SMIM26",
                        "HIST1H2AC" = "H2AC6",
                        "HIST1H2BC" = "H2BC4",
                        "DARC" = "ACKR1")
      data_combination <- data_combination %>%
        dplyr::mutate(GENE = recode(GENE, !!!replace_dict))
      df_name <- paste0(pheno, "_", dir, "_raw_value")
      assign(df_name, data_combination)
      print(paste0(df_name, " successfully generated"))
      if(keep_intermediate_file) {
        save_path <- file.path(intermediate_file_folder, paste0(df_name, ".txt"))
        write.table(data_combination, file=save_path, sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
      }
      ############### cs2g-step2:choose the highest linking in each sub strategy
      df <- data_combination
      #delete the non coding genes
      non_coding_genes <- c(
        "KRTAP5-AS1", "LOH12CR2", "LINC00332", "LINC00310", "SAMSN1-AS1",
        "RBM26-AS1", "LINC00160", "PPP1R2P2", "LINC00940", "CACNA1C-IT3",
        "RNASEH2B-AS1", "DLEU1", "GAS6-AS1", "LINC00652", "ADORA2A-AS1",
        "COL4A2-AS1", "LINC00163", "LINC00423", "MCF2L-AS1", "MCM3AP-AS1",
        "DIP2A-IT1", "LINC00598", "NPSR1-AS1", "NDFIP2-AS1", "RBMS3-AS3",
        "RPS16P5", "TRAM2-AS1", "EPHA1-AS1", "DLG5-AS1", "F10-AS1",
        "ZNRF3-AS1", "ARHGEF38-IT1", "MIR210HG", "KCNQ1-AS1", "DLEU2",
        "SUGT1P3", "ANKRD20A11P", "SLC22A18AS", "STARD13-AS", "SSR4P1",
        "CBR3-AS1", "TPTE2P5", "ITPK1-AS1", "LINC00567", "SMCR2",
        "CAPN15", "ALDH1L1-AS1", "MIR22HG", "ARHGEF3-AS1", "MAPT-AS1",
        "LINC00677", "LINC00605", "PDZRN3-AS1", "CLRN1-AS1", "SH3RF3-AS1",
        "ELFN2", "WEE2-AS1", "WWTR1-AS1", "TMEM105", "LINC00886",
        "TM4SF19-AS1", "TM4SF1-AS1", "POM121L9P", "BCRP3", "LRP5L",
        "RPL23AP82", "CRYBB2P1", "CABIN1", "UBXN7-AS1", "EZR-AS1",
        "LINC00624", "LINC01101", "PCGEM1", "HLTF-AS1", "LINC01003",
        "AIRN", "LINC00574", "MYLK-AS2", "C14orf178", "PSMB1",
        "LINC00303", "LINC01132", "LINC00184", "ALG1L9P", "SYNJ2-IT1",
        "TET2-AS1", "LL22NC01-81G9.3", "ZNF252P-AS1", "SPAG5-AS1", "SDCBP2-AS1",
        "PPP4R1L", "SMAD1-AS2", "SCARNA2", "SRRM2-AS1", "OVOL1-AS1",
        "CACTIN-AS1", "CASC8", "TNK2-AS1", "LINC00928", "KCNIP2-AS1",
        "VIM-AS1", "HOTAIRM1", "HOXA-AS2", "HOXA-AS3", "HOXA10-AS",
        "DBH-AS1", "LGALS8-AS1", "RPS2P32", "LINC00963", "EGFLAM-AS2",
        "SLC39A12-AS1", "LINC00682", "PRC1-AS1", "NR2F2-AS1", "LOXL1-AS1",
        "CASC19", "PCAT1", "CCDC26", "APOA1-AS", "LINC00861",
        "PVT1", "RNF138P1", "LINC00475", "CATSPER2P1", "HHIP-AS1",
        "ZNF252P", "CIRBP-AS1", "EGOT", "LIFR-AS1", "ADAMTSL4-AS1",
        "DIO3OS", "CASC11", "DLX6-AS1", "GAS5-AS1", "KIRREL3-AS3",
        "JAZF1-AS1", "KIAA0087", "CASC15", "IL21R-AS1", "RPS15AP10",
        "SNHG15", "TP53TG1", "UPK1A-AS1", "RUVBL1-AS1", "NCF1B",
        "UGDH-AS1", "DPY19L2P3", "LINC01094", "COLCA1", "P4HA2-AS1",
        "HPN-AS1", "LINC01063", "TIPARP-AS1", "LINC00486", "LINC00942",
        "ABCA17P", "FABP5P3", "ZNRF2P2", "CACNA1G-AS1", "ATP6V0E2-AS1",
        "SIGLEC17P", "APOC1P1", "DLG1-AS1", "ADCY10P1", "CCND2-AS1",
        "SCGB1B2P", "RAB11B-AS1", "LINC01001", "XXYLT1-AS1", "HNF1A-AS1",
        "ALDH1L1-AS2", "BHLHE40-AS1", "LINC00482", "GUCY2EP", "TPM3P9",
        "ZNF702P", "SEC1P", "MYLK-AS1", "PITPNA-AS1", "CSNK1G2-AS1",
        "LINC00877", "SEMA3B-AS1", "MBNL1-AS1", "LINC00881", "XXYLT1-AS2",
        "TMEM44-AS1", "SENCR", "SDHAP1", "RPL32P3", "FAM86JP",
        "FAM86HP", "PYY2", "ZNF286B", "LINC00674", "DNAH17-AS1",
        "LRP4-AS1", "MAPKAPK5-AS1", "LINC01023", "NCAM1-AS1", "KCNQ1OT1",
        "RAI1-AS1", "SMCR5", "NUP210P1", "DBIL5P", "SNHG16",
        "PCBP1-AS1", "SOX9-AS1", "MIR497HG", "FAHD2CP", "SMG7-AS1",
        "DNAJB3", "NEXN-AS1", "STARD7-AS1", "LINC00471", "TTN-AS1",
        "LINC00299", "PROX1-AS1", "RPL31P11", "PDIA3P1", "GNRHR2",
        "LEMD1-AS1", "SRP14-AS1", "CTD-3193O13.9", "FAM83H-AS1", "CSTF3-AS1",
        "CEBPA-AS1", "AATK-AS1", "GAS6-AS2", "TP73-AS1", "ITPR1-AS1", "SSSCA1-AS1",
        "TOLLIP-AS1", "KLHL7-AS1", "STT3A-AS1", "UBL7-AS1", "ANO1-AS2", "SPACA6P-AS",
        "UCHL1-AS1", "CDIPT-AS1", "BSN-AS2", "CHKB-AS1", "H1FX-AS1", "NCBP2-AS2",
        "PSMD5-AS1", "ACTN1-AS1", "TNRC6C-AS1", "C16orf80", "LINC01135", "C5orf56", "TMEM75", "C9orf139", "C17orf82", "CCAT1", "CRIPAK"
      )
      if(!is.null(gene_remove_list) && gene_remove_list == FALSE) {
        df <- df
      } else {
        genes_to_remove <- non_coding_genes
        if (!is.null(gene_remove_list)) {
          genes_to_remove <- unique(c(non_coding_genes, gene_remove_list))
        }
        df <- df[!(df$GENE %in% genes_to_remove), ]
      }

      result_df <- data.frame()

      for (i in 1:nrow(df)) {
        row <- df[i, ]
        if (is.na(row$VALUE)) {
          # if VALUE is NA, add directly
          result_df <- rbind(result_df, row)
        } else {
          # check whether there is any other rows as same as current SNP
          duplicate_rows <- df[df$SNP == row$SNP, ]
          if (nrow(duplicate_rows) == 1) {
            # if there are no duplicated SNP, set VALUE to 1
            row$VALUE <- 1
            result_df <- rbind(result_df, row)
          } else {
            # find the row with max VALUE in the duplicated SNP
            max_value <- max(duplicate_rows$VALUE, na.rm = TRUE)
            max_value_rows <- duplicate_rows[duplicate_rows$VALUE == max_value, ]

            if (row$VALUE == max_value) {
              # add the max VALUE column
              row$VALUE <- 1 / nrow(max_value_rows)
              result_df <- rbind(result_df, row)
            }
          }
        }
      }

      # check column name
      colnames(result_df) <- c("CNP", "GENE", "VALUE", "STRATEGY")
      result_df$SNP <- as.character(result_df$CNP)
      result_df$GENE <- as.character(result_df$GENE)
      result_df$VALUE <- as.numeric(result_df$VALUE)
      result_df$STRATEGY <- as.character(result_df$STRATEGY)
      # delete NA in VALUE column
      result_df <- result_df[!is.na(result_df$VALUE), ]
      # create object
      assign(paste(pheno, dir, "linking_score", sep = "_"), result_df, envir = .GlobalEnv)
      if(keep_intermediate_file) {
        save_path <- paste0(intermediate_file_folder, pheno, "_", dir, "_linking_score.txt")
        write.table(result_df, file=save_path, sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
      }
    }
    ##########cs2g step3:combine the link result from seven strategy
    if(!other_weights) {
      weights <- c(Exon=100, Promoter=100, finemappedciseQTLs_eQTLGen=25, finemappedciseQTLs_GTeX=7.5, EpiMap=1.9, ABC=1, Ciceroblood=1)
    }
    else {
      weights <- weights_vector
    }
    combined_df <- data.frame()
    for(dir in dir_list) {
      df_name <- paste(pheno, dir, "linking_score", sep = "_")
      df <- get(df_name)
      combined_df <- rbind(combined_df, df) %>% arrange(SNP)
    }
    ##Check if this pheno has one valid link at least
    if (sum(!is.na(combined_df$STRATEGY)) <= 0) {
      rlang::abort(paste("No successful linking in pheno:", pheno))
    }

    # Weighted Sum
    combined_df <- combined_df %>%
      group_by(SNP, GENE) %>%
      summarise(VALUE = sum(VALUE * sapply(STRATEGY, function(x) weights[x]), na.rm = TRUE),
                STRATEGY = toString(unique(STRATEGY)),
                .groups = 'drop')

    # rename
    combined_df <- combined_df %>%
      rename(SNP = SNP) %>%
      mutate(STRATEGY = stringr::str_replace_all(STRATEGY, "ciseQTLs_", ""))
    if(keep_intermediate_file) {
      save_path <- file.path(intermediate_file_folder,paste0(pheno,"_combined_score.txt"))
      write.table(combined_df,file=save_path,sep="\t",row.names=F,quote=F,col.names=T)
    }

    ########cs2g step4: standardized the linking score
    combined_df$VALUE_WEIGHTED <- NA
    # Iterate through each row to calculate VALUE_WEIGHTED
    for(i in 1:nrow(combined_df)) {
      current_snp <- combined_df$SNP[i]
      # Find all rows with the same SNP
      snp_rows <- combined_df[combined_df$SNP == current_snp, ]

      # Calculate the sum of the VALUE values for the duplicate rows
      sum_value <- sum(snp_rows$VALUE)

      # Update the VALUE_WEIGHTED column
      if(nrow(snp_rows) == 1) {
        # If there are no duplicate rows, set VALUE_WEIGHTED to 1
        combined_df$VALUE_WEIGHTED[i] <- 1
      } else {
        # If there are duplicate rows, calculate according to the formula
        combined_df$VALUE_WEIGHTED[i] <- round(combined_df$VALUE[i] / sum_value, 3)
      }
    }
    if(keep_intermediate_file) {
      save_path <- file.path(intermediate_file_folder,paste0(pheno,"_weighted_score.txt"))
      write.table(combined_df,file = save_path, sep="\t", row.names = F, quote = F, col.names = T)
    }
    ####### cs2g step5: calculate combined precision
    # filter the row with VALUE_WEIGHTED > 0.5
    filtered_df <- combined_df[combined_df$VALUE_WEIGHTED > 0.5, ]
    if (sum(!is.na(filtered_df$STRATEGY)) <= 0) {
      rlang::abort(paste("No successful linking in pheno:", pheno))
    }
    if(keep_intermediate_file) {
      save_path <- file.path(intermediate_file_folder,paste0(pheno,"_weighted_score_filtered.txt"))
      write.table(filtered_df,file = save_path, sep="\t", row.names = F, quote = F, col.names = T)
    }
    # precision map dict
    if(!other_precisions) {
      precision_values <- c(
        Exon = 1,
        finemappedGTeX = 0.676,
        EpiMap = 0.562,
        Promoter = 0.805,
        finemappedeQTLGen = 0.814,
        ABC = 0.469,
        Ciceroblood = 0.548
      )
    }
    else {
      precision_values <- other_precisions
    }

    calculate_combined_precision <- function(strategies) {
      precisions <- precision_values[strsplit(strategies, ",\\s*")[[1]]]
      combined_precision <- 1 - prod(1 - precisions)
      return(round(combined_precision, 3))
    }


    filtered_df$precision <- sapply(filtered_df$STRATEGY, calculate_combined_precision)
    if(causal_SNP) {
      filtered_df <- filtered_df %>%
        dplyr::left_join(SNP_data %>% select(SNP, PIP), by = "SNP") %>%
        mutate(confidence_score = precision * PIP)  # calculate 'confidence_score'
    }
    else {
      filtered_df$confidence_score <- filtered_df$precision * assumed_PIP
    }
    filtered_df <- filtered_df %>%
      dplyr::left_join(gene_annot_file, by = c("GENE" = "NAME"))
    if (!keep_unmatched_gene) {
      filtered_df <- filtered_df %>%
        filter(!is.na(ENSGID))
    }
    if (sum(!is.na(filtered_df$STRATEGY)) <= 0) {
      rlang::abort(paste("No successful linking after filtering genes that cannot be converted to ENSGID in pheno:", pheno))
    }
    write.csv(filtered_df,file.path(save_folder, paste0(pheno,"_cS2G.csv")),row.names=F)
  }

  if (length(pheno_name) > 1) {

    for (pheno in pheno_name) {
      current_SNP_data <- subset(SNP_data, SNP_data$PHENO == pheno)
      current_pheno_name <- unique(current_SNP_data$PHENO)
      tryCatch({
        link(
          pheno_name = current_pheno_name,
          SNP_data = current_SNP_data,
          s2g_folder = s2g_folder,
          save_folder = save_folder,
          keep_intermediate_file = keep_intermediate_file,
          intermediate_file_folder = intermediate_file_folder,
          gene_remove_list = gene_remove_list,
          other_weights = other_weights,
          other_precisions = other_precisions,
          causal_SNP = causal_SNP,
          assumed_PIP = assumed_PIP,
          gene_annot_file = gene_annot_file,
          keep_unmatched_gene = keep_unmatched_gene
        )
      }, error = function(e) {
        warning(paste("Skipping", pheno, "due to error:", e$message))
      })
    }
  } else {
    link(
      pheno_name = pheno_name,
      SNP_data = SNP_data,
      s2g_folder = s2g_folder,
      save_folder = save_folder,
      keep_intermediate_file = keep_intermediate_file,
      intermediate_file_folder = intermediate_file_folder,
      gene_remove_list = gene_remove_list,
      other_weights = other_weights,
      other_precisions = other_precisions,
      causal_SNP = causal_SNP,
      assumed_PIP = assumed_PIP,
      gene_annot_file = gene_annot_file,
      keep_unmatched_gene = keep_unmatched_gene
    )
  }
}




