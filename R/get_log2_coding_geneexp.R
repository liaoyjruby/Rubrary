#' Get log2p1 protein coding gene expression matrix
#'
#' @param exp_file Gene exp data; samples as columns, genes as rows
#' @param coding_file Protein coding genes mtx; "symbol" column
#' @param output_filename Filename to save output as
#'
#' @return txt (tsv) with log2p1 transformed protein coding only gene expression
#' @export
#'
get_log2_coding_geneexp <- function(exp_file, coding_file, output_filename){
  pc_genes <- read.delim(coding_file)
  pc_genes <- pc_genes$symbol
  exp <- read.delim(exp_file, row.names = 1)
  # Log2 transform
  exp_log2p1 <- log(exp + 1, base = 2)
  exp_log2p1_df <- tibble::rownames_to_column(as.data.frame(exp_log2p1), var = "gene")
  # Filter to coding only
  exp_log2p1_pc <- exp_log2p1_df[exp_log2p1_df$gene %in% pc_genes$symbol,]
  # Write to tsv
  write.table(exp_log2p1_pc,
              file = output_filename,
              sep = "\t", quote = F, row.names = F)
}
