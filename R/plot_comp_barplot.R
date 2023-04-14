#' Plot stacked barplot filled by group
#'
#' @import ggplot2
#'
#' @param sobj Seurat object
#' @param breaks string; metadata colname in sobj to set breaks / columns
#' @param group string; metadata colname in sobj to split columns by
#' @param stack logical; split into stacked plots of counts + fraction
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
plot_comp_barplot <- function(sobj, breaks, group, stack = TRUE,
                              title = ggplot2::waiver(), xlabel = breaks, colors = NULL,
                              ncol_legend = NULL, savename = NULL, width = 7, height = 7) {
  df <- sobj@meta.data[,c(breaks, group)]

  if(is.null(colors)){
    cols = scales::hue_pal()(length(unique(df[,group])))
  } else if (colors == "alpha"){
    cols = unname(pals::alphabet2(n = length(unique(df[,group]))))
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

    plt <- plt_ct / plt_scl +
      patchwork::plot_layout(heights = c(1,2))
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
