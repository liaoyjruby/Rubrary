#' Density plot by group
#'
#' @param df dataframe
#' @param value string; colname of values
#' @param group string; colname of group info
#' @param pval logical; get KSpval first group vs others
#' @param colors vector; # of colors = to # of groups
#' @param savename string; filepath to save PNG under
#'
#' @return Density plot by group
#' @export
#'
#' @examples
#' set.seed(13)
#'
plot_density <- function(df, value, group, pval = T,
                         colors = c("gray", "firebrick3"),
                         savename = NA) {
  # KS pval
  if (pval) {
    df <- df[order(df[, value]), ]
    df$rank <- 1:nrow(df)
    ks_pval <- stats::ks.test(
      df[df[, group] == unique(df[, group])[1], "rank"],
      df[!(df[, group] == unique(df[, group])[1]), "rank"]
    )$p.value
  }

  dplt <- ggpubr::ggdensity(
    data = df,
    x = value,
    xlab = "",
    ylab = "Density",
    fill = group,
    alpha = 0.85,
    palette = colors
  ) +
    theme_classic() + {
      if (pval) labs(subtitle = paste0("KS enrich. p-value = ", signif(ks_pval, digits = 4)))
    }

  if (!is.na(savename)) {
    ggsave(
      filename = savename,
      plot = dplt,
      width = 6, height = 4
    )
  }

  return(dplt)
}
