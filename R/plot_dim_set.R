
#' Plot Seurat DimPlot grid set
#'
#' @param sobj Seurat object
#' @param reduction c("pca", "tsne", "umap); must be already run in sobj
#' @param set dataframe or vector; col1 metadata name, col2 (opt) subplot description, col3 (opt) logical T/F for label
#' @param ncol number of columns in grid
#' @param guides string; specify how guides should be treated in layout (?patchwork::wrap_plots)
#' @param redtitle logical; `TRUE` to prepend reduction method in subplot title
#' @param savename string; filepath to save figure under
#' @param width integer; width of plot
#' @param height integer; height of plot
#'
#' @return Patchwork grid of Seurat DimPlots
#' @export
plot_dim_set <- function(sobj, reduction = "umap", set, ncol = 3,
                         guides = NULL, redtitle = FALSE,
                         savename = NULL, width = NULL, height = NULL){
  plot_dim <- function(group, name, lb){
    red <- ifelse(reduction == "tsne", "tSNE", toupper(reduction))
    if(methods::is(sobj@meta.data[,group], "numeric")){
      plt <- Seurat::FeaturePlot(
        sobj,
        reduction = reduction,
        features = group,
        cols = c("blue", "red")
      ) +
        ggplot2::labs(title = ifelse(redtitle, paste0(red, " - ", name), name))
    } else {
      plt <- Seurat::DimPlot(
        sobj,
        reduction = reduction,
        group.by = group,
        shuffle = T,
        label = lb,
        label.box = lb
      ) +
        ggplot2::labs(title = ifelse(redtitle, paste0(red, " - ", name), name))
    }
    return(plt)
  }

  if(is.vector(set)){
    set <- data.frame(group = set)
  }
  if(ncol(set) == 1) {set[,2] = set[,1]}
  if(ncol(set) == 2){ set[,3] = FALSE}

  plts <- mapply(plot_dim, set[,1], set[,2], set[,3], SIMPLIFY = FALSE)
  dim_set <- patchwork::wrap_plots(plts, ncol = ncol, guides = guides)

  if(!is.null(savename)) {
    ggplot2::ggsave(
      filename = savename,
      plot = dim_set,
      width = ifelse(is.null(width), ncol * 7 + 2, width),
      height = ifelse(is.null(height), ceiling(nrow(set) / ncol) * 6, height)
    )
  }
  return(dim_set)
}
