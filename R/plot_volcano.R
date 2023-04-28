#' Plot volcano plot
#'
#' Wrapper for `EnhancedVolcano::EnhancedVolcano()` with additional caption and saving
#'
#' @param df_deg dataframe; test statistics, w/ cols "names", "x", "y"
#' @param names string; colname in `df_deg` w/ variable names
#' @param x string; colname for log2 fold changes
#' @param y string; colname for nominal/adj p-values
#' @param pCutoff numeric; value to set p value cutoff
#' @param FCcutoff numeric; value to set log fold change cutoff
#' @param title string; plot title
#' @param subtitle string; plot subtitle
#' @param xlab_high string; description of high x values
#' @param xlab_low string; description of low x values
#' @param xlab_size numeric; size of x labels
#' @param savename string; full name of file to save under
#' @param height numeric; saved plot height
#' @param width numeric; saved plot width
#'
#' @return ggplot2 volcano plot
#' @export
#'
plot_volcano <- function(df_deg, names = "gene", x = "log2FoldChange", y = "pvalue",
                         pCutoff = 5e-4, FCcutoff = 1,
                         xlab_high = NULL, xlab_low = NULL, xlab_size = 7,
                         title = NULL, subtitle = NULL,
                         savename = NULL, height = 12, width = 16){
  Rubrary::use_pkg("EnhancedVolcano")

  ypos_top = round(max(-log10(df_deg[,y])) * 0.8, 2)

  xlims <- c(min(df_deg[[x]], na.rm = TRUE) - 1.5,
             max(df_deg[[x]], na.rm = TRUE) + 1.5)
  xpos_left = round(xlims[1] + (0.5 * abs(xlims[1])), 2)
  xpos_right = round(xlims[2] - (0.5 * xlims[2]), 2)

  vol_plt <- EnhancedVolcano::EnhancedVolcano(
    df_deg,
    lab = df_deg[,names],
    x = x,
    y = y,
    FCcutoff = FCcutoff,
    pCutoff = pCutoff,
    labSize = 4,
    pointSize = 1.5,
    title = title,
    subtitle = subtitle,
    caption = paste0(y, " cutoff = ", pCutoff, "; ", x," cutoff = ", FCcutoff)
  ) +
    {if (!is.null(xlab_high)) ggplot2::geom_text(
      x = xpos_right, y = ypos_top, label = xlab_high, size = xlab_size, hjust = 0.5)} +
    {if (!is.null(xlab_low)) ggplot2::geom_text(
      x = xpos_left, y = ypos_top, label = xlab_low, size = xlab_size, hjust = 0.5)}

  if (!is.null(savename)){
    ggplot2::ggsave(
      filename = savename,
      plot = vol_plt,
      height = 12,
      width = 16
    )
  }

  return(vol_plt)
}
