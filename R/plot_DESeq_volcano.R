#' Plot DESeq volcano
#'
#' @param df_deseq dataframe; DESeq output read as df
#' @param pCutoff numeric; value to set p value cutoff
#' @param logFCcutoff numeric; value to set log fold change cutoff
#' @param title string; plot title
#' @param subtitle string; plot subtitle
#' @param savename string; full name of file to save under
#'
#' @return DE genes volcano plot
#' @export
#'
plot_DESeq_volcano <- function(df_deseq,
                               pCutoff = 5e-4, logFCcutoff = 1,
                               title = "DE Genes", subtitle = "", savename = NA){

  vol_plt <- EnhancedVolcano::EnhancedVolcano(
    df_deseq,
    lab = df_deseq$gene,
    x = 'log2FoldChange',
    y = 'pvalue',
    FCcutoff = logFCcutoff,
    pCutoff = pCutoff,
    labSize = 4,
    pointSize = 1.5,
    title = title,
    subtitle = subtitle,
    caption = paste0("p-value cutoff = ", pCutoff, "; logFC cutoff = ", logFCcutoff)
  )

  if (!is.na(savename)){
    ggplot2::ggsave(
      filename = savename,
      plot = vol_plt,
      height = 12,
      width = 16
    )
  }

  return(vol_plt)
}
