#' plot_boxplot_paired
#'
#' @param sig Column name of specific signature
#' @param sig_paired_df Dataframe with all signatures
#' @param batch Plot subtitle
#'
#' @return Paired boxplot separated by type
#' @export
#'
plot_boxplot_paired <- function(sig, sig_paired_df, batch){
  sig_paired_df <- sig_paired_df[,c("Sample","Type", sig)]
  # order <- c(rbind(1:(nrow(sig_paired_df)/2), ((nrow(sig_paired_df)/2)+1):nrow(sig_paired_df)))
  # sig_paired_df <- sig_paired_df[order,] # Put pairs next to each other
  colnames(sig_paired_df) <- c("Sample", "Type", "Sum.Z.Score")
  ggpubr::ggpar(
    ggpubr::ggpaired(
      sig_paired_df,
      x = "Type",
      y = "Sum.Z.Score",
      xlab = F,
      ylab = "Sum Z Score",
      title = sig,
      subtitle = batch,
      color = "black",
      fill = "Type",
      line.color = "gray",
      line.size = 0.4,
      palette = "jco",
      label = sig_paired_df$Sample,
      font.label = list(size = 8, color = "black"),
      repel = T
    ),
    legend = "none"
  ) +
    ggpubr::stat_compare_means(paired = TRUE, label = "p.format", label.x = 1.3)
}
