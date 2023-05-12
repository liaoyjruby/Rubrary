#' Plot simple scatter
#'
#' Correlation value and method gets placed in plot caption / bottom-right.
#'
#' @import ggplot2
#'
#' @param df dataframe; includes both signatures
#' @param xval string/numeric vector; colname for x axis
#' @param yval string/numeric vector; colname for y axis
#' @param label string; colname for point labels
#' @param group string; colname for group color
#' @param title string; plot title
#' @param subtitle string; plot subtitle
#' @param rank logical; change metric to rank
#' @param xlabel string; xval description
#' @param ylabel string; yval description
#' @param cormethod string; correlation method for stats::cor() or "none"
#' @param pt_size numeric; point size
#' @param pt_alpha numeric; point alpha value
#' @param lbl_size numeric; label text size
#' @param guides logical; T to include fitted linear model line
#' @param reverse logical; T to reverse both axes
#' @param heatmap logical; T for underlying heatmap (untested)
#' @param hm_palette string; RColorBrewer continuous palette name (untested)
#' @param savename string; filepath to save figure under
#' @param width numeric; plot width
#' @param height numeric; plot height
#'
#' @return Simple scatter plot as ggplot2
#' @export
#'
#' @examples
#' df <- data.frame(A = 1:10, B = 11:20, C = LETTERS[1:10], D = rep(LETTERS[1:5], each = 2))
#' plot_scatter(df, xval = "A", yval = "B", label = "C", group = "D", title = "Example Scatter")
#'
plot_scatter <- function(df = NULL, xval, yval, label = NULL, group = NULL, rank = FALSE,
                         xlabel = ifelse(methods::is(xval, "character"), xval, "X"),
                         ylabel = ifelse(methods::is(yval, "character"), yval, "Y"),
                         cormethod = c("pearson", "spearman", "none"), guides = TRUE,
                         pt_size = 2, pt_alpha = 1, lbl_size = 3,
                         heatmap = FALSE, hm_palette = "Spectral",
                         reverse = FALSE, title = NULL, subtitle = NULL,
                         savename = NULL, width = 8, height = 8) {
  cormethod <- match.arg(cormethod)

  if(is.null(df)){
    if (!is.null(label)){
      df <- data.frame(xval = xval, yval = yval, label = label)
    } else {
      df <- data.frame(xval = xval, yval = yval)
    }
    xval <- "xval"
    yval <- "yval"
    label <- "label"
  } else {
      df <- df[, c(xval, yval, label, group)]
  }

  if (rank) {
    df <- df[order(df[,xval], decreasing = T),]
    df[,xval] <- 1:nrow(df)
    df <- df[order(df[,yval], decreasing = T), ]
    df[,yval] <- 1:nrow(df)
    cormethod <- "spearman"
    # corcoef <- "rho"
    corcoef <- "\U03C1"
    guides <- FALSE
    xlabel <- paste0(xlabel, " Rank")
    ylabel <- paste0(ylabel, " Rank")
  } else {
    corcoef <- "R"
  }

  if(cormethod != "none"){
    corr <- signif(cor(
      x = df[, xval],
      y = df[, yval],
      method = cormethod
    ), 2)
  }


  plt <- ggplot(df, aes(x = .data[[xval]], y = .data[[yval]])) +
    {if(guides) geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red")} +
    {if(heatmap) geom_raster()} +
    {if(heatmap) scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, hm_palette)),
      breaks = scales::pretty_breaks(5)) } +
    geom_point(alpha = pt_alpha, size = pt_size) +
    {if(!is.null(group)) geom_point(aes(color = .data[[group]]), alpha = pt_alpha, size = pt_size)} +
    xlab(xlabel) +
    ylab(ylabel) +
    labs(
      title = title,
      subtitle = subtitle
    ) +
    {if(cormethod != "none") labs(caption = paste0(
      corcoef, " = ", corr, "; ", tools::toTitleCase(cormethod), " correlation"))} +
    theme_classic() +
    {if(reverse) scale_x_reverse()} +
    {if(reverse) scale_y_reverse()} +
    {if(!is.null(label)) ggrepel::geom_text_repel(aes(label = .data[[label]]), max.overlaps = Inf, size = lbl_size)}

  if (!is.null(savename)) {
    ggplot2::ggsave(
      filename = savename,
      plot = plt,
      height = height,
      width = width
    )
  }

  return(plt)
}
