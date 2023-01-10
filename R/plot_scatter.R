#' Plot simple scatter
#'
#' @param df Dataframe; includes both signatures
#' @param xval String; colname for x axis
#' @param yval String; colname for y axis
#' @param eqlims Equal limits? TO BE IMPLEMENTED
#' @param save Logical; save as .png
#' @param outname Output file name
#'
#' @return Simple scatter plot
#' @importFrom ggplot2 ggplot aes_string geom_point geom_smooth theme_classic
#' @export
#'
plot_scatter <- function(df = NA, xval, yval, eqlims = F, save = F, outname = "Scatter.png") {
  df <- df[,c(xval, yval)]

  if(eqlims){
    xrange <- c(min(df[,xval], na.rm = T), max(df[,xval], na.rm = T))
    yrange <- c(min(df[,yval], na.rm = T), max(df[,yval], na.rm = T))
    limits <- c(min(xrange, yrange) * 1.1, max(xrange, yrange) * 1.1)
  }

  plt <- ggplot(df, aes_string(x = xval, y = yval)) +
    geom_point() +
    geom_smooth(method = "lm", se=FALSE, color = 'red') +
    ggpubr::stat_cor(method = "pearson", label.x.npc = 0.4) +
    theme_classic()

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
