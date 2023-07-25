
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

#' Plot stacked barplot filled by group
#'
#' @import ggplot2
#'
#' @param sobj Seurat object
#' @param breaks string; metadata colname in sobj to set breaks / columns
#' @param group string; metadata colname in sobj to split columns by
#' @param stack logical; split into stacked plots of counts + fraction
#' @param counts logical; T to include counts barplot above stacked barplot
#' @param title string; plot title
#' @param xlabel string; x-axis (breaks) label
#' @param colors char vector; list of colors (n = # of unique groups), or "alpha"
#' @param ncol_legend integer; # of legend columns
#' @param savename string; filepath to save plot under
#' @param width numeric; saved plot width
#' @param height numeric; saved plot height
#'
#' @return Stacked barplot, filled by group composition
#' @export
#'
plot_comp_barplot <- function(sobj, breaks, group, stack = TRUE, counts = TRUE,
                              title = ggplot2::waiver(), xlabel = breaks, colors = NULL,
                              ncol_legend = NULL, savename = NULL, width = 7, height = 7) {
  df <- sobj@meta.data[,c(breaks, group)]

  if(is.null(colors)){
    Rubrary::use_pkg("scales")
    cols = scales::hue_pal()(length(unique(df[,group])))
  } else if (colors[1] == "alpha"){
    Rubrary::use_pkg("Seurat")
    cols = Seurat::DiscretePalette(n = length(unique(df[,group])), palette = "alphabet2")
  } else {
    cols = colors
  }

  if (!stack){
    # Unscaled
    plt <- ggplot(df, aes(x = .data[[breaks]], fill = .data[[group]])) +
      geom_bar(stat="count") +
      labs(title = title) +
      ylab("# of cells") +
      xlab(xlabel) +
      scale_fill_manual(values=cols) +
      theme_classic() +
      guides(fill = guide_legend(ncol=ncol_legend)) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      ) + patchwork::plot_layout() # Helps scale legend appropriately
  } else {
    # Split and scaled
    plt_ct <- ggplot(df, aes(x = .data[[breaks]])) +
      geom_bar(position = "stack", stat="count", fill = "black") +
      labs(title = title) +
      ylab("# of cells") +
      theme_classic() +
      guides(fill = guide_legend(ncol=ncol_legend)) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank()
      )

    plt_scl <- ggplot(df, aes(x = .data[[breaks]], fill = .data[[group]])) +
      geom_bar(position = "fill", stat="count") +
      ylab("% of cells") +
      xlab(xlabel) +
      scale_fill_manual(values=cols) +
      theme_classic() +
      guides(fill = guide_legend(ncol=ncol_legend)) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )

    if(counts){
      plt <- plt_ct / plt_scl +
        patchwork::plot_layout(heights = c(1,2))
    } else {
      plt <- plt_scl
    }
  }

  if(!is.null(savename)){
    ggsave(
      filename = savename,
      plot = plt,
      width = width, height = height
    )
  }
  return(plt)
}

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
  corr_method <- match.arg(corr_method)

  Rubrary::use_pkg("Seurat")
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
