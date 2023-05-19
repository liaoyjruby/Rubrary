utils::globalVariables(c(
  "cms", "cms_smooth"
))

#' Compute and plot integration assessment cell-specific methods for Seurat object
#'
#' Wrapper for application of single cell integration assessment metrics with default parameters. Takes in Seurat object, calculates integration assessment metric per cell for given grouping variables, and plots visualization of metrics. A single integration score per dataset can be reported as the median of the metric across all cells.
#'
#' ## Local inverse Simpson's index `lisi`: neighborhood batch diversity
#'
#' Wrapper for [immunogenomics/LISI](https://github.com/immunogenomics/LISI) from Harmony paper. Quantifies diversity of batches within local neighborhoods, defined by kNN graphs, using the inverse Simpson's index.
#' - LISI = 1: neighborhood comprised of single group, not mixed
#' - LISI = # of groups: neighborhood has equal number of cells from all datasets, perfect mixing
#'
#' ## Seurat MixingMetric `mm`: local / group kNN ranking
#'
#' Wrapper for `Seurat::MixingMetric` from Seurat v3 paper (Stuart et al. 2019). Assesses local neighborhood composition for a cell and how well mixed that neighborhood is. If the neighborhood is well mixed, should contain at least some cells from each dataset. Can be more intuitive to have high MixingMetric mean good mixing, so range can be flipped with `max.k - MixingMetric`. Uses flipped scale (`MM_flip = TRUE`) by default.
#' - MixingMetric = `max.k`, low `max.k - MixingMetric`: not mixed
#' - Low MixingMetric, high `max.k - MixingMetric`: well mixed
#' - Set-wide metric: median (or mean) of MixingMetric distribution
#'
#' Implemented with `max.k = 300` by default.
#'
#' ## CellMixS `cms`: batch-specific distance distribution kNN test
#'
#' Wrapper for `CellMixS::cms` (Lütge et al. 2021). Tests for each cell the hypothesis that batch-specific distance distributions towards it's k-nearest neighboring (knn) cells are derived from the same unspecified underlying distribution using the Anderson-Darling test.
#' - Low `cms`: poorly mixed per cell neighborhood
#' - High `cms`: well mixed per cell neighborhood
#' - Set-wide metric: `cms` histogram, flat for random batch mixing and skewed towards 0 for batch related bias
#'
#' Implementation uses `k = 300` and `dim_red = reduction`. As per [the CellMixS vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/CellMixS/inst/doc/CellMixS.html#parameter), `k` should not exceed the size of the smallest cell population. Takes a while to run, so asks for `BiocParallel` to use `BiocParallel::MulticoreParam()`. Computer may churn while it runs...
#'
#' @import dplyr
#'
#' @param sobj Seurat object
#' @param group.var character (vector for `mm` or `lisi`); metadata colname to check level of integration
#' @param method c("LISI", "MixingMetric", "CellMixS"); integration assessment metric to apply
#' @param reduction c("umap", "pca", "tsne"); dimensional reduction method
#' @param title string; plot title
#' @param savename string; filepath to save results and figure under (no ext.)
#' @param k integer; `max.k` for `Seurat::MixingMetric`, `k` for `CellMixS::cms`
#' @param width numeric; plot width
#' @param height numeric; plot height
#' @param MM_flip logical; use `max.k - MixingMetric` instead of `MixingMetric` only
#'
#' @return Dataframe of integration assessment values
#' @seealso [lisi::compute_lisi()], [Seurat::MixingMetric()], [CellMixS::cms()]
#' @export
#'
assess_integration <- function(sobj, group.var, method = c("LISI", "MixingMetric", "CellMixS"),
                     reduction = c("umap", "pca", "tsne"), k = 300, MM_flip = TRUE,
                     title = NULL, savename = NULL, width = NULL, height = NULL){
  Rubrary::use_pkg("Seurat")
  reduction <- match.arg(reduction)
  method <- match.arg(method)
  wd <- width
  ht <- height

  # Resolve missing values
  meta <- sobj@meta.data %>%
    select(all_of(group.var)) %>%
    mutate(across(everything(), ~ tidyr::replace_na(.x, "NA")))

  if(method == "MixingMetric"){ # Seurat::MixingMetric
    max.k = k # Default params

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
      xlabel <- "Seurat MixingMetric (mm, flipped)"
      interpret <- "Seurat MixingMetric (mm, flipped) interpretation:\n** Low `max.k - mm`: poorly mixed per cell neighborhood\n** High `max.k - mm`: well mixed per cell neighborhood\n** Set-wide metric: median (or mean) of `max.k - mm` distribution"
    } else { # Use MixingMetric
      colnames(res) <- paste0("MixingMetric_", group.var)
      xlabel <- "Seurat MixingMetric (mm)"
      interpret <- "Seurat MixingMetric (mm) interpretation:\n** High `mm`: poorly mixed per cell neighborhood\n** Low `mm`: well mixed per cell neighborhood\n** Set-wide metric: median (or mean) of `mm` distribution"
    }
  } else if (method == "CellMixS") {
    Rubrary::use_pkg("SingleCellExperiment", "CellMixS", "BiocParallel")

    # Convert to SingleCellExperiment
    sce <- Seurat::as.SingleCellExperiment(sobj)

    ndims <- ncol(SingleCellExperiment::reducedDim(sce, toupper(reduction)))

    sce <- CellMixS::cms(
      sce = sce,
      k = k,
      group = group.var,
      dim_red = toupper(reduction),
      n_dim = ndims,
      BPPARAM = BiocParallel::MulticoreParam()
    )

    res <- SingleCellExperiment::colData(sce) %>% # metadata df
      as.data.frame() %>%
      select(cms, cms_smooth)
    colnames(res) <- c(paste0("CellMixS_", group.var), paste0("CellMixS_smooth_", group.var))

    xlabel <- "Cell-specific Mixing Score (cms)"
    interpret <- "CellMixS (cms) interpretation:\n** Low `cms` = 0: poorly mixed per cell neighborhood\n** High `cms` = 1: well mixed per cell neighborhood\n** Set-wide metric: `cms` histogram, flat for random batch mixing and skewed towards 0 for batch related bias"

  } else { # LISI by default
    Rubrary::use_pkg("lisi")
    res <- lisi::compute_lisi(
      X = as.data.frame(Seurat::Embeddings(sobj, reduction = reduction)),
      meta_data = meta,
      label_colnames = group.var)
    colnames(res) <- paste0("LISI_", group.var)

    xlabel <- "Local Inverse Simpson's Index (lisi)"
    interpret <- "Local Inverse Simpson's Index (lisi) interpretation:\n** Min `lisi` = 1: cell neighborhood consists of one group, poorly mixed\n** Max `lisi` = number of groups: cell neighborhood consists of all groups, well mixed\n** Set-wide metric: median (or mean) of `lisi` distribution"
  }

  message(interpret) # interpretation blurb

  if (!is.null(savename)){
    utils::write.table(
      x = tibble::rownames_to_column(res, var = "Sample"),
      file = paste0(savename, "_results.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }

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
      wd <- 16
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
        hist = (method == "CellMixS")
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

  if(!is.null(savename)){
    ggplot2::ggsave(
      plot = grid,
      filename = paste0(savename, "_grid.png"),
      width = wd, height = ht
    )
  }

  return(res)
}

#' @describeIn assess_integration Compute and plot local inverse Simpson's index `lisi`
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

#' @describeIn assess_integration Compute and plot Seurat MixingMetric `mm`
#' @export
run_MixingMetric <- function(sobj, group.var, MM_flip = TRUE, reduction = c("umap", "pca", "tsne"),
                             max.k = 300, title = NULL, savename = NULL, width = NULL, height = NULL){
  Rubrary::assess_integration(
    sobj = sobj, method = "MixingMetric",
    MM_flip = MM_flip,
    k = max.k,
    group.var = group.var,
    reduction = reduction,
    title = title,
    savename = savename,
    width = width, height = height
  )
}

#' @describeIn assess_integration Compute and plot CellMixS `cms`
#' @export
run_CellMixS <- function(sobj, group.var, reduction = c("umap", "pca", "tsne"),
                         k = 300, title = NULL, savename = NULL, width = NULL, height = NULL){
  Rubrary::assess_integration(
    sobj = sobj, method = "CellMixS",
    k = k,
    group.var = group.var,
    reduction = reduction,
    title = title,
    savename = savename,
    width = width, height = height
  )
}
