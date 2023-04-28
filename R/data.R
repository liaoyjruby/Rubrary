#' Airway Dex Treated vs. Untreated DESeq2 Results
#'
#' `airway` package data run through DESeq2 to analyze differential gene expression between dexamethasone treated and untreated airway smooth muscle cells. Original gene names converted from Ensembl ID to HGNC symbol and filtered to protein-coding only through the BioMart database.
#'
#' @format ## `airway_deseq_res`
#' A data frame with 16,460 rows and 9 columns:
#' \describe{
#'   \item{ensembl_gene_id, hgnc_symbol}{Gene names in Ensembl ID and HGNC symbol format}
#'   \item{baseMean}{Average of normalized count values, dividing by size factors, taken over all samples}
#'   \item{log2FoldChange, lfcSE}{Effect size estimate, on log base 2 scale and its standard error}
#'   \item{stat}{Value of test statistic}
#'   \item{pvalue, padj}{P-value and adjusted P-value}
#'   \item{sign_log_p}{Log base 2 of p-value multiplied by sign}
#'   ...
#' }
#' @source <https://bioconductor.org/packages/release/data/experiment/html/airway.html>
"airway_deseq_res"
