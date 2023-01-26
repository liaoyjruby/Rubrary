#' Output DESeq results as text and rank file
#'
#' @param DE_results DESeq results object
#' @param outname string; output filename sans extenstion
#'
#' @return DESeq results as text and rank file
#' @export
#'
output_DESeq <- function(DE_results, outname = "DESeq_results") {
  res_df <- as.data.frame(DE_results)
  res_df <- res_df[!(res_df$baseMean == 0),] # Filter out all where mean was 0
  res_df$sign_log_p <- sign(res_df$stat)*-log2(res_df$pvalue) # Calc signed log p
  res_df <- res_df[!is.na(res_df$sign_log_p),] # NA slogp unimportant
  res_df <- res_df[order(res_df$sign_log_p, decreasing = T),] # Reorder
  res_out <- tibble::rownames_to_column(res_df, "gene")
  write.table(x = res_out,
              file = paste0(outname, ".txt"),
              sep = "\t", row.names = F, quote = F)
  # Rank output
  res_rank <- res_out[,c("gene", "sign_log_p")]
  res_rank$sign_log_p_rank <- 1:nrow(res_rank)
  write.table(res_rank[,c("gene", "sign_log_p_rank")],
              file = paste0(outname, ".rnk"),
              sep = '\t', quote = F, row.names = F)

  return(res_out)
}
