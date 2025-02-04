% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eqtl_coloc.R
\name{eqtl_coloc}
\alias{eqtl_coloc}
\title{Colocalization-based linking}
\usage{
eqtl_coloc(
  pheno_name,
  SNP_data_create_above,
  gwas_separate_data_folder,
  eqtl_separate_data_folder,
  strand_harmonize,
  gwas_data_type,
  gwas_sdy_input,
  gwas_sdy_data = NULL,
  gwas_s_data = NULL,
  gwas_beta_input,
  gwas_n_data = NULL,
  eqtl_n_data = NULL,
  save_folder,
  work_part = NULL,
  coloc_p1 = NULL,
  coloc_p2 = NULL,
  coloc_p12 = NULL
)
}
\arguments{
\item{coloc_p12}{}
}
\description{
This function utilizes "coloc::coloc.abf" to perform colocalization analysis between eQTL and GWAS, providing colocalization-based linking evidence.
}
