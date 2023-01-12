#' Plot simple scatter
#'
#' @param df Dataframe; includes both signatures
#' @param xval String; colname for x axis
#' @param yval String; colname for y axis
#' @param label String; colname for point labels
#' @param eqlims Equal limits? TO BE IMPLEMENTED
#' @param save Logical; save as .png
#' @param outname Output file name
#' @param title String; title of plot
#'
#' @return Simple scatter plot
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth theme_classic labs
#' @export
#'
#' @examples
#' df <- data.frame(A = 1:10, B = 11:20, C = LETTERS[1:10])
#' plot_scatter(df, xval = "A", yval = "B", label = "C", title = "Example Scatter")
#'
plot_scatter <- function(df = NA, xval, yval, label = NA, eqlims = F,
                         title = paste0(xval, " vs. ", yval),
                         save = F, outname = "Scatter.png") {
  df <- df[, c(xval, yval, label)]

  if (eqlims) {
    xrange <- c(min(df[, xval], na.rm = T), max(df[, xval], na.rm = T))
    yrange <- c(min(df[, yval], na.rm = T), max(df[, yval], na.rm = T))
    limits <- c(min(xrange, yrange) * 1.1, max(xrange, yrange) * 1.1)
  }

  plt <- ggplot(df, aes_string(x = xval, y = yval, label = label)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    ggpubr::stat_cor(method = "pearson", label.x.npc = 0.4) +
    labs(
      title = title,
    ) +
    theme_classic() +
    {
      if (!is.na(label)) ggrepel::geom_text_repel(max.overlaps = Inf)
    }

  if (save) {
    ggplot2::ggsave(
      filename = outname,
      plot = plt,
      height = 8,
      width = 8
    )
  }

  return(plt)
}
