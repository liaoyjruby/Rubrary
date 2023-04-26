utils::globalVariables(c(
  "baseMean", "pvalue", "sign_log_p"
))

#' Output DESeq results as dataframe with signed log2 p metric
#'
#' Filters out all `sign_log_p == NA` rows.
#'
#' @import dplyr
#'
#' @param DE_results DESeq results object
#' @param savename string; filepath to save table under (no ext.)
#' @param rank logical; T to output .rnk file compatible with Java GSEA app
#'
#' @return Dataframe with DESeq results
#' @export
#'
output_DESeq <- function(DE_results, savename = NULL, rank = FALSE) {
  res_df <- as.data.frame(DE_results) %>%
    filter(baseMean != 0) %>% # Filter out all where mean was 0
    mutate(sign_log_p = sign(stat)*-log2(pvalue)) %>% # Calc signed log p
    filter(!is.na(sign_log_p)) %>%
    arrange(desc(sign_log_p)) %>%
    tibble::rownames_to_column(var = "gene")

  if(!is.null(savename)){
    write.table(x = res_df,
                file = paste0(savename, ".txt"),
                sep = "\t", row.names = F, quote = F)
  }

  if(rank){
    res_rank <- res_df %>%
      select(c(gene, sign_log_p)) %>%
      mutate(sign_log_p_rank = 1:nrow(res_df))
    if(!is.null(savename)){
      write.table(res_rank[,c("gene", "sign_log_p_rank")],
                  file = paste0(savename, ".rnk"),
                  sep = '\t', quote = F, row.names = F)
    } else {
      res_df <- res_rank
    }
  }
  return(res_df)
}
