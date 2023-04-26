#' Plot simple scatter
#'
#' @import ggplot2
#'
#' @param df dataframe; includes both signatures
#' @param xval string/numeric vector; colname for x axis
#' @param yval string/numeric vector; colname for y axis
#' @param label string; colname for point labels
#' @param savename string; filepath to save figure under
#' @param title string; plot title
#' @param rank logical; change metric to rank
#' @param xlabel string; xval description
#' @param ylabel string; yval description
#' @param cormethod string; correlation method for stats::cor()
#' @param alpha numeric; point alpha value
#' @param guides logical; T to include fitted line
#' @param reverse logical; T to reverse both axes
#' @param heatmap logical; T for underlying heatmap (untested)
#' @param hm_palette string; RColorBrewer continuous palette name (untested)
#'
#' @return Simple scatter plot as ggplot2
#' @export
#'
#' @examples
#' df <- data.frame(A = 1:10, B = 11:20, C = LETTERS[1:10])
#' plot_scatter(df, xval = "A", yval = "B", label = "C", title = "Example Scatter")
#'
plot_scatter <- function(df = NULL, xval, yval, label = NULL, rank = FALSE,
                         xlabel = ifelse(methods::is(xval, "character"), xval, "X"),
                         ylabel = ifelse(methods::is(yval, "character"), xval, "Y"),
                         cormethod = c("pearson", "spearman"), guides = TRUE,
                         alpha = 1, heatmap = FALSE, hm_palette = "Spectral",
                         reverse = FALSE, title = paste0(xval, " vs. ", yval),
                         savename = NULL) {
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
    if (!is.null(label)){
      df <- df[, c(xval, yval, label)]
    } else {
      df <- df[, c(xval, yval)]
    }
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

  # if (eqlims) {
  #   xrange <- c(min(df[, xval], na.rm = T), max(df[, xval], na.rm = T))
  #   yrange <- c(min(df[, yval], na.rm = T), max(df[, yval], na.rm = T))
  #   limits <- c(min(xrange, yrange) * 1.1, max(xrange, yrange) * 1.1)
  # }

  corr <- signif(cor(
    x = df[, xval],
    y = df[, yval],
    method = cormethod
  ), 2)

  plt <- ggplot(df, aes(x = .data[[xval]], y = .data[[yval]])) +
    {if(guides) geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red")} +
    {if(heatmap) geom_raster()} +
    {if(heatmap) scale_fill_gradientn(
      colors = rev(RColorBrewer::brewer.pal(11, hm_palette)),
      breaks = scales::pretty_breaks(5)) } +
    geom_point(alpha = alpha, size = 0.5) +
    # ggpubr::stat_cor(method = cormethod, cor.coef.name = corcoef) +
    xlab(xlabel) +
    ylab(ylabel) +
    labs(
      title = title,
      subtitle = paste0(corcoef, " = ", corr, "; ", tools::toTitleCase(cormethod), " correlation")
    ) +
    theme_classic() +
    {if(reverse) scale_x_reverse()} +
    {if(reverse) scale_y_reverse()} +
    {if(!is.null(label)) ggrepel::geom_text_repel(aes(label = .data[[label]]), max.overlaps = Inf)}

  if (!is.null(savename)) {
    ggplot2::ggsave(
      filename = savename,
      plot = plt,
      height = 8,
      width = 8
    )
  }

  return(plt)
}
