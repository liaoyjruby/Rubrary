#' Plot simple scatter
#'
#' @param df Dataframe
#' @param xval Colname for x axis
#' @param yval Colname for y axis
#' @param eqlims Equal limits?
#'
#' @return Simple scatter plot
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth
#' @export
#'
plot_scatter <- function(df, xval, yval, eqlims = F) {
  df <- df[,c(xval, yval)]

  if(eqlims){
    xrange <- c(min(df[,xval], na.rm = T), max(df[,xval], na.rm = T))
    yrange <- c(min(df[,yval], na.rm = T), max(df[,yval], na.rm = T))
    limits <- c(min(xrange, yrange) * 1.1, max(xrange, yrange) * 1.1)
  }

  ggplot(df, aes_string(x = xval, y = yval)) +
    geom_point() +
    geom_smooth(method = "lm", se=FALSE, color = 'red') +
    ggpubr::stat_cor(method = "pearson", label.x.npc = 0.4)
}
