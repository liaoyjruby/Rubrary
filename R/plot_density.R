#' Density plot by group
#'
#' @param df dataframe
#' @param value string; colname of values
#' @param group string; colname of group info
#' @param group2 string; group of interest
#' @param pval logical; get KSpval first group vs others
#' @param colors vector; # of colors = to # of groups
#' @param savename string; filepath to save PNG under
#' @param title string; plot title
#' @param xlab string; x-axis label
#' @param rug logical; put rug plot on x axis
#' @param rug_label logical; label rug ticks with first df column
#' @param rug_lab_nudge numeric; fine-tune distance of label from x-axis
#'
#' @return Density plot by group
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_rug
#' @export
#'
#' @examples
#' set.seed(13)
#'
plot_density <- function(df, value, group, group2 = NA, title = NA, pval = T,
                         xlab = value, colors = c("firebrick3", "gray"),
                         rug = F, rug_label = F, rug_lab_nudge = 0.00001,
                         savename = NA,width = 6, height = 4) {
  if(is.na(group2)) {
    group2 <- unique(df[,group])[1]
  }

  # KS pval
  if (pval) {
    ks_pval <- Rubrary::get_kspval(df, value, group, group2)
  }

  dplt <- ggplot(data = df,
                 aes(x = .data[[value]], fill = .data[[group]])) +
    {if(rug) geom_rug(data = df[df[,group] == group2,], color = colors[1])} +
    geom_density(alpha = 0.85) +
    scale_fill_manual(values = colors) +
    {if(rug_label) ggrepel::geom_text_repel(
      data = df[df[,group] == group2,],
      mapping = aes(x = .data[[value]], y = 0, label = df[df[,group] == group2,][,1],),
      angle = 90, hjust = 0, direction = 'x', nudge_y = rug_lab_nudge)} +
    ylab("Density") +
    xlab(xlab) +
    {if (!is.na(title)) labs(title = title)} +
    {if (pval) labs(subtitle = paste0("KS enrich. p-value = ", signif(ks_pval, digits = 4)))} +
    theme_classic()

  if (!is.na(savename)) {
    ggsave(
      filename = savename,
      plot = dplt,
      width = width, height = height
    )
  }

  return(dplt)
}
