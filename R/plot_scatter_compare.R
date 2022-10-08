
#' plot_scatter_compare
#'
#' @param set1path Path to 1st set (x axis)
#' @param set2path Path to 2nd set (y axis)
#' @param set1lab Axis label for 1st set (x)
#' @param set2lab Axis label for 2nd set (y)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param rank Convert sign log p values to rank
#' @param outname File output
#' @param save Save as .png
#'
#' @return Plot
#' @export
#'
plot_scatter_compare <- function(set1path, set2path, set1lab = "Set 1", set2lab = "Set 2",
                                 rank = F, save = F,
                                 title = paste0(set1lab, " vs. ", set2lab), subtitle = "",
                                 outname = paste0(dirname(set1path), "/", set1lab, "vs", set2lab, ".png")) {

  sig1 <- utils::read.delim(set1path, header = T)
  sig2 <- utils::read.delim(set2path, header = T)

  if (rank) {
    sig1 <- sig1[order(sig1$sign_log_p, decreasing = T), ]
    sig2 <- sig2[order(sig2$sign_log_p, decreasing = T), ]
    sig1$sign_log_p <- 1:nrow(sig1)
    sig2$sign_log_p <- 1:nrow(sig2)
  }
  sig_merged <- base::merge(sig1[, c("gene", "sign_log_p")], sig2[, c("gene", "sign_log_p")], by = "gene", all = T)
  colnames(sig_merged) <- c("gene", "sign_log_p.1", "sign_log_p.2")
  sig1_range <- c(min(sig_merged$sign_log_p.1, na.rm = T), max(sig_merged$sign_log_p.1, na.rm = T))
  sig2_range <- c(min(sig_merged$sign_log_p.2, na.rm = T), max(sig_merged$sign_log_p.2, na.rm = T))
  limits <- c(min(sig1_range, sig2_range) * 1.1, max(sig1_range, sig2_range) * 1.1)
  plt <- ggplot2::ggplot(sig_merged, ggplot2::aes_string(x = "sign_log_p.1", y = "sign_log_p.2")) +
    ggplot2::geom_point(alpha = 0.2, size = 0.5) +
    # geom_abline(linetype="dashed", ggplot2::aes(intercept=0, slope=1), size = 1) +
    # geom_smooth(method = "lm", se=FALSE, color = 'red') +
    ggpubr::stat_cor(
      method = "pearson",
      label.x = limits[1] + (limits[2] * 0.05),
      label.y = limits[2]
    ) +
    # stat_regline_equation(label.x = -40, label.y = 28) +
    # coord_fixed(xlim = limits, ylim = limits) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title, # include comp % on axis
      subtitle = subtitle,
      x = set1lab,
      y = set2lab
    )

  if (save) {
    ggplot2::ggsave(
      filename = outname,
      plot = plt,
      height = 8,
      width = 8
    )
  }

  plt
}
