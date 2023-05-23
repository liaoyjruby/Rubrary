
#' Plot Seurat DimPlot grid set
#'
#' @import ggplot2
#'
#' @param sobj Seurat object
#' @param reduction c("pca", "tsne", "umap); must be already present in sobj
#' @param set dataframe or vector; col1 metadata name, col2 (opt) subplot description, col3 (opt) logical T/F for label
#' @param cols_num char vector; two colors for low numeric legend range to high numeric legend range
#' @param ncol number of columns in grid
#' @param guides string; specify how guides should be treated in layout (?patchwork::wrap_plots)
#' @param redtitle logical; `TRUE` to prepend reduction method in subplot title
#' @param gridtitle string; title for overall dimplot grid
#' @param gridsubtitle string; subtitle for overall dimplot grid
#' @param savename string; filepath to save figure under
#' @param width integer; plot width, if unspecified will try to calculate
#' @param height integer; plot height, if unspecified will try to calculate
#'
#' @return Patchwork grid of Seurat DimPlots
#' @export
plot_dimplot_grid <- function(sobj, reduction = c("umap", "pca", "tsne"), set, ncol = 3,
                         guides = NULL, redtitle = FALSE, gridtitle = NULL, gridsubtitle = NULL,
                         cols_num = c("blue", "red"),
                         savename = NULL, width = NULL, height = NULL){
  Rubrary::use_pkg("Seurat")

  reduction <- match.arg(reduction)

  plot_dim <- function(group, name, lb){
    red <- ifelse(reduction == "tsne", "tSNE", toupper(reduction))
    if(!(group %in% names(sobj@meta.data)) ||
       (methods::is(sobj@meta.data[,group], "numeric"))){
      plt <- Seurat::FeaturePlot(
        sobj,
        reduction = reduction,
        features = group,
        cols = cols_num
      ) +
        labs(title = ifelse(redtitle, paste0(red, " - ", name), name))
    } else {
      plt <- Seurat::DimPlot(
        sobj,
        reduction = reduction,
        group.by = group,
        shuffle = T,
        label = lb,
        label.box = lb
      ) +
        labs(title = ifelse(redtitle, paste0(red, " - ", name), name))
    }
    return(plt)
  }

  if(is.vector(set)){
    set <- data.frame(group = set)
  }
  if(ncol(set) == 1) {set[,2] = set[,1]}
  if(ncol(set) == 2){ set[,3] = FALSE}

  plts <- mapply(plot_dim, set[,1], set[,2], set[,3], SIMPLIFY = FALSE)
  dim_set <- patchwork::wrap_plots(plts, ncol = ncol, guides = guides) +
    patchwork::plot_annotation(
      title = gridtitle,
      subtitle = gridsubtitle)

  if(!is.null(savename)) {
    ggsave(
      filename = savename,
      plot = dim_set,
      width = ifelse(is.null(width), ncol * 7 + 2, width),
      height = ifelse(is.null(height), ceiling(nrow(set) / ncol) * 6, height)
    )
  }
  return(dim_set)
}
