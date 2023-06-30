#' Get log(x+1) protein coding gene expression matrix
#'
#' @import dplyr
#' @param exp df/string; samples as columns, genes as rows
#' @param base integer; log transformation base
#' @param pc_genes char vector; list of coding genes to filter by
#' @param savename string; filename to save output as
#'
#' @return Log(x+1) transformed protein coding only gene expression dataframe
#' @export
#'
get_log_coding_geneexp <- function(exp, base = 2, pc_genes = Rubrary::get_PC_genes(),
                                   savename){
  if(is.character(exp)){ exp <- Rubrary::rread(exp, row.names = 1) }
  # Log transform
  exp_logp1_pc <- log(exp + 1, base = base) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "gene") %>%
    Rubrary::filter_genes(genes_filt = pc_genes, gene_col = "gene")
  # Write to tsv
  Rubrary::rwrite(
    x = exp_logp1_pc,
    file = savename
  )
  return(exp_logp1_pc)
}
