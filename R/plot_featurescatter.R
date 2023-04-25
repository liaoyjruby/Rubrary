#' Plot FeatureScatter wrapper
#'
#' Wrapper for `Seurat::FeatureScatter` with a few extras like saving options, plot title, and correlation in the plot subtitle
#'
#' @param sobj Seurat object
#' @param feat1 string; 1st feature to plot (x-axis)
#' @param feat2 string; 2nd feature to plot (y-axis)
#' @param group string; metadata colname to group by
#' @param title string; plot title
#' @param savename string; filepath to save plot under
#' @param corr_method string; `stats::cor` method for correlation calculation
#' @param height numeric; plot height
#' @param width numeric; plot width
#'
#' @return FeatureScatter ggplot output
#' @export
plot_featurescatter <- function(sobj, feat1, feat2, group = "sample_id", title = NULL,
                                savename = NULL, corr_method = c("pearson", "kendall", "spearman"),
                                height = 7, width = 9) {
  Rubrary::use_pkg("Seurat")
  corr_method <- match.arg(corr_method)

  corr <- round(stats::cor(x = Seurat::FetchData(sobj, feat1),
                    y = Seurat::FetchData(sobj, feat2),
                    method = corr_method), digits = 2)

  p <- Seurat::FeatureScatter(
    object = sobj,
    feature1 = feat1,
    feature2 = feat2,
    shuffle = T,
    group.by = group
  ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0(tools::toTitleCase(corr_method), " correlation = ", corr)
    ) +
    ggplot2::theme_classic()

  if(!is.null(savename)){
    ggplot2::ggsave(
      plot = p,
      filename = savename,
      width = width, height = height
    )
  }
  return(p)
}
