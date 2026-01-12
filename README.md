# mdS2G: Multi-Dimensional SNP-to-Gene Linking Framework

## Overview

**mdS2G** is a R package designed to prioritize causal genes by integrating multiple lines of evidence:

1.  **Locus-based linking** (based on cS2G)
2.  **Similarity-based linking** (PoPS)
3.  **Colocalization-based linking** (Coloc)

This tool was developed for the paper *"Cross-ancestry GWAS of liver function-related biomarkers enhance insights into causal variants, credible genes, and therapeutic targets"*. It provides a robust pipeline for linking GWAS signals to their credible causal genes across different phenotypes. This package is  designed to run in a **Linux environment**.

**Note:** Some external databases (e.g., cS2G annotations, eQTL summary statistics) are too large to be included in the package. Download links for these required resources are provided in the relevant sections below.

## System Requirements

To ensure smooth execution, your system should meet the following requirements:

* **R Version**: >= 4.3
* **Package Dependencies**:
    * AnnotationDbi (>= 1.64.0)
    * data.table (>= 1.17.6)
    * dplyr (>= 1.1.4)
    * magrittr
    * org.Hs.eg.db (>= 3.18.0)
    * R.utils
    * rlang
    * stats
    * stringr
    * tidyr (>= 1.3.1)
    * utils

## Installation

You can install the development version of `mdS2G` from GitHub using the `devtools` package.

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("Y-Wentao/mdS2G")

# Create a directory for testing
if (!dir.exists("mds2g_test")) {
  dir.create("mds2g_test")
}
setwd("mds2g_test")

# Load the library
library(mds2g)

```

## Step-by-Step Usage Guide

This demo uses **ALT** and **AST** (liver function traits) as example phenotypes to demonstrate the workflow.

### Step 1: Locus-based Linking (cS2G)

Perform fine-mapping based annotation using the modified cS2G framework.

**Input Data Requirements**: Prior to running this function, you must conduct *in-silico* fine-mapping and retain putative causal SNPs with a Posterior Inclusion Probability (**PIP**) $\ge$ 0.5. The input data frame must contain the columns `CHR`, `SNP`, `BP`, and `PIP`. When analyzing multiple phenotypes simultaneously, a `PHENO` column is also required to identify the phenotype corresponding to each SNP.

For further details, please consult the associated manuscript and the original cS2G publication (DOI: 10.1038/s41588-022-01087-y).

**Prerequisite**: Download the cS2G annotation file (`S2G_bed.tgz`) from [Zenodo](https://zenodo.org/records/10117202) and unzip it.

```r
cs2g_link(
  pheno_name = c("ALT", "AST"),
  SNP_data = "./input/fine_mapping_results.csv",  # Path to your fine-mapping file
  s2g_folder = "/path/to/downloaded/S2G_bed/",    # Path to the unzipped S2G_bed folder
  save_folder = "./output/",                      # Directory to save results
  keep_intermediate_file = TRUE,
  intermediate_file_folder = "./output/"
)

```

**Output**: `ALT_cS2G.csv` and `AST_cS2G.csv` in the output folder.

### Step 2: Similarity-based Linking (PoPS) Processing

Process results from PoPS (Polygenic Priority Score).

**Prerequisite**: Before executing this step, you must independently perform PoPS analysis using the same GWAS results to generate the final `.preds` files (e.g., `ALT.preds`, `AST.preds`). As PoPS annotations rely on Ensembl Gene IDs (ENSGID), this function is specifically tailored to handle this identifier format. For implementation details, please refer to the original PoPS study (DOI: 10.1038/s41588-023-01443-6).

```r
pops_result_processing(
  pheno_name = c("ALT", "AST"),
  SNP_data = "./input/fine_mapping_results.csv",  # Must match input from Step 1
  pops_result_folder = "./input/pops_results/",   # Folder containing your .preds files
  save_folder = "./output/"
)

```

**Output**: `*_pops_result_processed.csv` files.

### Step 3: Nominated Gene Extraction

Integrate results from Step 1 & 2 to nominate genes for colocalization and filter high-confidence genes.

```r
nominated_gene_extract(
  pheno_name = c("ALT", "AST"),
  cs2g_result_folder = "./output/",   # Result folder from Step 1
  pops_result_folder = "./output/",   # Result folder from Step 2
  save_folder = "./output/"
)

```

**Output**:

* `gene_for_coloc.txt`: List of genes for eQTL extraction.
* `SNP_data_for_coloc.csv`: Loci list for colocalization.
* `cs2g_cred_gene.csv`: High-confidence credible genes from cS2G.

### Step 4: eQTL Data Preparation

Clean and format cis-eQTL data (SMR format) for colocalization. 
**Note:** You may choose to use data from all available tissues or select a specific subset relevant to your study. This function is specifically optimized for the data sources listed below due to their high accessibility and standardized formatting. **Crucially, the input files must be kept in their original compressed state (e.g., .tar.gz) without any modification or decompression.**
**If you intend to use alternative data versions, other tissues, or different types of QTL results for colocalization, you must manually format your data to mirror the output structure of this function (i.e., QTL summary statistics must be split by both tissue and gene).**

**Prerequisite**:

1. **Data**: Download SMR-formatted cis-eQTL data for [GTEx v8](https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary.html) or [eQTLGen](https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz).
2. **Software**: Download [SMR executable](https://yanglab.westlake.edu.cn/software/smr/#Download).

```r
eqtl_data_preparation(
  eqtl_data_folder = "/path/to/downloaded/eQTL_data/", # Folder with downloaded eQTL data
  smr_exe = "/path/to/software/smr",                   # Path to SMR executable file
  extract_gene_list = "./output/gene_for_coloc.txt"    # Generated in Step 3
)

```

**Output**: Separate folders for each tissue containing split eQTL data.

### Step 5: GWAS Data Preparation

Format GWAS summary statistics for colocalization.
**Input Data Requirements**: The input GWAS summary statistics **must** contain `SNP`, `CHR`, and `BP` columns. Additional columns are required depending on the specific parameters selected. For instance, in this example, the columns `BETA`, `SE`, `N`, `MAF`, `A1`, and `A2` are explicitly required. **Please ensure that the column names in your input files match these names exactly.**

```r
gwas_data_preparation(
  pheno_name = c("ALT", "AST"),
  SNP_data_create_above = "./output/SNP_data_for_coloc.csv", # Generated in Step 3
  gwas_data_folder = "./input/gwas_summary_stats/",          # Folder with your GWAS sumstats
  gwas_data_suffix = ".txt.gz",
  gwas_allele_input = TRUE,
  gwas_data_type = "quant",       # "quant" or "cc"
  gwas_sdy_input = FALSE,         # Set TRUE if sdY is available
  gwas_beta_input = TRUE,         # Set TRUE if beta is available
  save_folder = "./output/gwas_prepared/"
)

```

**Output**: Folders for each phenotype containing locus-specific summary statistics.

### Step 6: Colocalization Analysis

Run `coloc.abf()` to test for colocalization between GWAS and eQTL signals.

```r
eqtl_coloc(
  pheno_name = c("ALT", "AST"),
  SNP_data_create_above = "./output/SNP_data_for_coloc.csv",
  gwas_separate_data_folder = "./output/gwas_prepared/",
  eqtl_separate_data_folder = "./output/eqtl_prepared/", # Output from Step 4
  strand_harmonize = TRUE,
  gwas_data_type = "quant",
  gwas_sdy_input = FALSE,
  gwas_beta_input = TRUE,
  gwas_n_data = "./input/gwas_samplesize.csv", # Required for quant traits without sdY.
  save_folder = "./output/coloc_results/"
)

```

**Output**:

* `coloc.csv`: All locus-gene pairs with colocalization evidence.
* `coloc_SNP.csv`: Specific colocalized SNPs.

### Step 7: Final Result Processing

Integrate all evidence to report the final list of **Credible Genes**.

```r
result_process(
  coloc_result = "./output/coloc_results/coloc.csv",
  SNP_data_create_above = "./output/SNP_data_for_coloc.csv",
  cs2g_cred_data = "./output/cs2g_cred_gene.csv",
  save_folder = "./output/final_results/"
)

```

**Final Output Files**:

1. `coloc_result_processed.csv`: Cleaned colocalization results.
2. `credible_gene_muti_evidence.csv`: Genes supported by multiple strategies (e.g., PoPS + Coloc).
3. `CredibleGene_all.csv`: Comprehensive list of all identified credible genes.
4. `unique_CredibleGene.csv`: A matrix of unique credible genes across phenotypes.

## References

If you use **mdS2G** in your research, please cite our manuscript:

* *Cross-ancestry GWAS of liver function-related biomarkers enhance insights into causal variants, credible genes, and therapeutic targets.* (Under Review)

**Please also cite the original publications for the underlying methodologies if they were incorporated into your analysis:**

* **cS2G**: [DOI: 10.1038/s41588-022-01087-y](https://doi.org/10.1038/s41588-022-01087-y)
* **PoPS**: [DOI: 10.1038/s41588-023-01443-6](https://doi.org/10.1038/s41588-023-01443-6)
* **Coloc**: [DOI: 10.1371/journal.pgen.1004383](https://doi.org/10.1371/journal.pgen.1004383)
* **SMR**: [DOI: 10.1038/ng.3538](https://doi.org/10.1038/ng.3538)

## Contact

For any questions or issues, please contact:
wentao.yao@foxmail.com
