
#' Plot scatter, either metric or rank
#'
#' Rubrary::plot_scatter is written to be a bit more general.
#'
#' @import ggplot2
#'
#' @param set1path Path to 1st set (x axis)
#' @param set2path Path to 2nd set (y axis)
#' @param set1lab Axis label for 1st set (x)
#' @param set2lab Axis label for 2nd set (y)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param rank Convert sign log p values to rank
#' @param savename string; File output full name
#' @param guides logical; Show m=1 and linear fit?
#'
#' @return Plot comparing DESeq signatures
#' @export
#'
plot_scatter_compare <- function(set1path, set2path, set1lab = "Set 1", set2lab = "Set 2", rank = F, guides = F,
                                 title = paste0(set1lab, " vs. ", set2lab), subtitle = "",
                                 savename = NULL) {

  sig1 <- utils::read.delim(set1path, header = T)
  sig2 <- utils::read.delim(set2path, header = T)

  sig_merged <- base::merge(sig1[, c("gene", "sign_log_p")], sig2[, c("gene", "sign_log_p")], by = "gene")
  colnames(sig_merged) <- c("gene", "sign_log_p.1", "sign_log_p.2")
  if (rank) {
    sig_merged <- sig_merged[order(sig_merged$sign_log_p.1, decreasing = T), ]
    sig_merged$sign_log_p.1 <- 1:nrow(sig_merged)
    sig_merged <- sig_merged[order(sig_merged$sign_log_p.2, decreasing = T), ]
    sig_merged$sign_log_p.2 <- 1:nrow(sig_merged)
    cormethod = "spearman"
  } else {
    cormethod = "pearson"
  }
  sig1_range <- c(min(sig_merged$sign_log_p.1, na.rm = T), max(sig_merged$sign_log_p.1, na.rm = T))
  sig2_range <- c(min(sig_merged$sign_log_p.2, na.rm = T), max(sig_merged$sign_log_p.2, na.rm = T))
  limits <- c(min(sig1_range, sig2_range) * 1.1, max(sig1_range, sig2_range) * 1.1)
  plt <- ggplot(sig_merged, aes_string(x = "sign_log_p.1", y = "sign_log_p.2")) +
    geom_point(alpha = 0.2, size = 0.5) +
    {if (!rank) geom_abline(linetype="dashed", aes(intercept=0, slope=1), size = 1)} +
    {if (!rank) geom_smooth(method = "lm", se=FALSE, color = 'red')} +
    ggpubr::stat_cor(method = cormethod) +
    theme_classic() +
    labs(
      title = title, # include comp % on axis
      subtitle = paste0(subtitle, "; ", tools::toTitleCase(cormethod), " correlation"),
      x = set1lab,
      y = set2lab
    )

  if (!is.null(savename)) {
    ggsave(
      filename = savename,
      plot = plt,
      height = 8,
      width = 8
    )
  }

  return(plt)
}
