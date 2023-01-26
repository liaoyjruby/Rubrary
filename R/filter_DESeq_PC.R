#' Filter DESeq signature to protein coding only
#'
#' @param sig_txt Path to DESeq text file output
#' @param coding_file Path to protein-coding gene table w/ "symbol" column
#'
#' @return DESeq text and rnk files subsetted to PC coding only gene
#' @importFrom utils read.delim
#' @importFrom tools file_path_sans_ext
#' @export
#'
filter_DESeq_PC <- function(sig_txt, coding_file = "/Users/liaoyj/Dropbox/Ovarian Project/log2_coding_expression_datasets/protein-coding_gene.txt") {
  pc_genes <- read.delim(coding_file, header = T)
  sig_full <- read.delim(sig_txt, header = T)
  sig_full_pc <- sig_full[sig_full$gene %in% pc_genes$symbol, ]
  write.table(x = sig_full_pc,
              file = paste0(file_path_sans_ext(sig_txt), "_PC.txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  # Make rank
  sig_rnk_pc <- sig_full_pc[,c("gene", "sign_log_p")]
  write.table(sig_rnk_pc,
              file = paste0(file_path_sans_ext(sig_txt), "_PC.rnk"),
              sep = "\t", quote = F, col.names = F, row.names = F)

  return(sig_full_pc)
}
