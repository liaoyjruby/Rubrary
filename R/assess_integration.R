#' Compute and plot integration assessment methods for Seurat object
#'
#' Wrapper for application of single cell integration assessment metrics with default parameters. Takes in Seurat object, calculates integration assessment metric per cell for given grouping variables, and plots visualization of metrics. A single integration score per dataset can be reported as the median of the metric across all cells.
#'
#' ## LISI: local inverse Simpson's index
#'
#' Wrapper for [immunogenomics/LISI](https://github.com/immunogenomics/LISI) from Harmony paper. Quantifies diversity of batches within local neighborhoods, defined by kNN graphs, using the inverse Simpson's index.
#' - LISI = 1: neighborhood comprised of single group, not mixed
#' - LISI = # of groups: neighborhood has equal number of cells from all datasets, perfect mixing
#'
#' ## Seurat MixingMetric: local / group kNN ranking
#'
#' Wrapper for `Seurat::MixingMetric` from Seurat v3 paper. Assesses local neighborhood composition for a cell and how well mixed that neighborhood is. If the neighborhood is well mixed, should contain at least some cells from each dataset. Can be more intuitive to have high MixingMetric mean good mixing, so range can be flipped with `max.k - MixingMetric`. Uses flipped scale (`MM_flip = TRUE`) by default.
#' - MixingMetric = `max.k`, low `max.k - MixingMetric`: not mixed
#' - Low MixingMetric, high `max.k - MixingMetric`: well mixed
#'
#'
#' @import dplyr
#'
#' @param sobj Seurat object
#' @param group.var character (vector); metadata colname to check level of integration
#' @param method c("LISI", "MixingMetric"); integration assessment metric to apply
#' @param reduction c("umap", "pca", "tsne"); dimensional reduction method
#' @param title string; plot title
#' @param savename string; filepath to save results and figure under (no ext.)
#' @param width numeric; plot width
#' @param height numeric; plot height
#' @param MM_flip logical; use `max.k - MixingMetric` instead of `MixingMetric` only
#'
#' @return Dataframe of integration assessment values
#' @seealso [Seurat::MixingMetric()]
#' @export
#'
assess_integration <- function(sobj, group.var, method = c("LISI", "MixingMetric"),
                     reduction = c("umap", "pca", "tsne"), MM_flip = TRUE,
                     title = NULL, savename = NULL, width = NULL, height = NULL){
  Rubrary::use_pkg("Seurat", "lisi")

  reduction <- match.arg(reduction)
  method <- match.arg(method)
  wd <- width
  ht <- height

  # Resolve missing values
  meta <- sobj@meta.data %>%
    select(all_of(group.var)) %>%
    mutate(across(everything(), ~ tidyr::replace_na(.x, "NA")))

  if(method == "MixingMetric"){ # Seurat::MixingMetric
    max.k = 300 # Default params

    res <- lapply(
      group.var,
      function(gvar){
        Seurat::MixingMetric(
          object = sobj,
          grouping.var = gvar,
          reduction = reduction,
          max.k = max.k)}) %>% as.data.frame()
    rownames(res) <- colnames(sobj)
    if(MM_flip){ # Use max.k - MixingMetric
      colnames(res) <- paste0("MixingMetric_", group.var, "_orig")
      mm_rev <- max.k - res[,paste0("MixingMetric_", group.var, "_orig")] %>%
        as.data.frame()
      colnames(mm_rev) <- paste0("MixingMetric_", group.var)
      res <- cbind(mm_rev, res)
      xlabel <- "Seurat MixingMetric (flipped)"
    } else {
      colnames(res) <- paste0("MixingMetric_", group.var)
      xlabel <- "Seurat MixingMetric"
    }
  } else { # LISI by default
    res <- lisi::compute_lisi(
      X = as.data.frame(Seurat::Embeddings(sobj, reduction = reduction)),
      meta_data = meta,
      label_colnames = group.var)
    colnames(res) <- paste0("LISI_", group.var)

    xlabel <- "Local Inverse Simpson's Index (LISI)"
  }

  if (!is.null(savename)){
    utils::write.table(
      x = tibble::rownames_to_column(res, var = "Sample"),
      file = paste0(savename, "_results.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }

  if(!is.null(savename)){
    # Add metric to metadata
    sobj <- Seurat::AddMetaData(object = sobj, metadata = res)

    # Smart (?) grid calcs
    if(length(group.var) > 1){
      nc <- length(group.var)
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

    # Metric dimplot set
    dimset <- data.frame(
      name = c(group.var, paste0(method, "_", group.var)),
      desc = c(group.var, paste0(toupper(reduction), " - ", group.var, " ", method)),
      label = rep(FALSE, length(c(group.var, paste0(method, "_", group.var))))
    )

    dimplots <- Rubrary::plot_dimplot_grid(
      sobj = sobj,
      reduction = reduction,
      set = dimset,
      ncol = nc
    )

    # Metric distribution density plot
    dists <- lapply(
      paste0(method, "_", group.var),
      function(x){
        Rubrary::plot_distribution(
          values = res[,x],
          title = paste0(method, " Distribution"),
          xlab = xlabel,
        )
      }
    )
    distplots <- patchwork::wrap_plots(dists, ncol = length(dists))

    if(length(group.var) >1){
      grid <- (dimplots / distplots) + patchwork::plot_layout(heights = c(8,1)) +
        patchwork::plot_annotation(title = title)
    } else {
      grid <- (dimplots | distplots) + patchwork::plot_layout(widths = c(2,1)) +
        patchwork::plot_annotation(title = title)
    }

    print(grid)

    ggplot2::ggsave(
      plot = grid,
      filename = paste0(savename, "_grid.png"),
      width = wd, height = ht
    )
  }

  return(res)
}

#' @describeIn assess_integration Compute and plot local inverse Simpson's index (LISI)
#' @export
run_LISI <- function(sobj, group.var, reduction = c("umap", "pca", "tsne"),
                     title = NULL, savename = NULL, width = NULL, height = NULL){
  Rubrary::assess_integration(
    sobj = sobj, method = "LISI",
    group.var = group.var,
    reduction = reduction,
    title = title,
    savename = savename,
    width = width, height = height
  )
}

#' @describeIn assess_integration Compute and plot Seurat MixingMetric
#' @export
run_MixingMetric <- function(sobj, group.var, MM_flip = TRUE, reduction = c("umap", "pca", "tsne"),
                     title = NULL, savename = NULL, width = NULL, height = NULL){
  Rubrary::assess_integration(
    sobj = sobj, method = "MixingMetric",
    MM_flip = MM_flip,
    group.var = group.var,
    reduction = reduction,
    title = title,
    savename = savename,
    width = width, height = height
  )
}
