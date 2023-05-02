#' Compute and plot local inverse Simpson index
#'
#' Wrapper for [immunogenomics/LISI](https://github.com/immunogenomics/LISI) that takes in Seurat object, calculates LISI per cell for given labels/variables, and plots visualization of LISI results.
#'
#' @import dplyr
#'
#' @param sobj Seurat object
#' @param labels character vector; variable colname in metadata to calculate LISI values against
#' @param reduction c("umap", "pca", "tsne"); dimensional reduction method
#' @param title string; plot title
#' @param savename string; filepath to save results and figure under (no ext.)
#' @param width numeric; plot width
#' @param height numeric; plot height
#'
#' @return Dataframe of LISI values
#' @export
#'
run_LISI <- function(sobj, labels, reduction = c("umap", "pca", "tsne"),
                     title = NULL, savename = NULL, width = NULL, height = NULL){
  Rubrary::use_pkg("Seurat", "lisi")

  reduction <- match.arg(reduction)
  wd <- width
  ht <- height

  # Resolve missing values
  meta <- sobj@meta.data %>%
    select(all_of(labels)) %>%
    mutate(across(everything(), ~ tidyr::replace_na(.x, "NA")))

  # Run LISI
  res <- lisi::compute_lisi(
    X = as.data.frame(Seurat::Embeddings(sobj, reduction = reduction)),
    meta_data = meta,
    label_colnames = labels)
  colnames(res) <- paste0("lisi_", labels)

  # Add LISI to metadata
  sobj <- Seurat::AddMetaData(object = sobj, metadata = res)

  if(length(labels) > 1){
    nc <- length(labels)
    if(is.null(wd) || is.null(ht)){
      wd <- nc * 6 + 2
      ht <- 15
    }
  } else { # Single row plot
    nc <- 2
    if(is.null(wd) || is.null(ht)){
      wd <- 15
      ht <- 5
    }
  }

  # LISI dimplot set
  lisi_set <- data.frame(
    name = c(labels, paste0("lisi_", labels)),
    desc = c(labels, paste0(
      labels, " miLISI = ",
      unlist(lapply(paste0("lisi_", labels),
                    function(x){round(stats::median(res[,x]), digits = 2)})))),
    label = rep(FALSE, length(c(labels, paste0("lisi_", labels))))
  )

  lisi_dimplots <- Rubrary::plot_dimplot_grid(
    sobj = sobj,
    reduction = reduction,
    set = lisi_set,
    ncol = nc
  )

  # LISI distribution
  lisi_dists <- lapply(
    paste0("lisi_", labels),
    function(x){
      Rubrary::plot_distribution(
        values = res[,x],
        title = x,
        xlab = "Local Inverse Simpson's Index (LISI)",
      )
    }
  )
  lisi_distplots <- patchwork::wrap_plots(lisi_dists, ncol = length(lisi_dists))

  if(length(labels) >1){
    lisi_grid <- (lisi_dimplots / lisi_distplots) + patchwork::plot_layout(heights = c(8,1)) +
      patchwork::plot_annotation(title = title)
  } else {
    lisi_grid <- (lisi_dimplots | lisi_distplots) + patchwork::plot_layout(widths = c(2,1)) +
      patchwork::plot_annotation(title = title)
  }

  if(!is.null(savename)){
    utils::write.table(
      x = tibble::rownames_to_column(res, var = "Sample"),
      file = paste0(savename, "_results.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    ggplot2::ggsave(
      plot = lisi_grid,
      filename = paste0(savename, "_grid.png"),
      width = wd, height = ht
    )
  }

  print(lisi_grid)

  return(res)
}
