#' Plot GSEA pathway plots (waterfall + enrichment)
#'
#' @import ggplot2
#'
#' @param sig dataframe; "gene" + rankcol columns
#' @param rankcol string; colname of values
#' @param rankcol_name string; descriptive name of values
#' @param geneset vector; list of genes
#' @param title string; title of plot
#' @param subtitle string; (name of signature?)
#' @param savename string; filepath to save png under
#' @param highlab string; description of high values
#' @param lowlab string; description of low values
#' @param hllab string; description of highlighted values
#'
#' @return grid of gene set enrichment waterfall above enrichment plot
#' @export
#'
plot_GSEA_pathway <- function(sig, geneset, rankcol, rankcol_name = rankcol,
                              highlab = NA, lowlab = NA, hllab = "Highlight",
                              title = "", subtitle = NA, savename = NULL){
  path_genes <- geneset
  label <- length(path_genes) < 20
  # Waterfall plot
  plt_wf <- Rubrary::plot_waterfall(
    sig = sig,
    label = label,
    highlight = path_genes,
    highlab = highlab,
    lowlab = lowlab,
    hllab = hllab, # hllab = "Highlight",
    otherlab = "Other genes", # otherlab = "Other",
    rankcol = rankcol,
    ylab = rankcol_name,
    title = title
  ) +
    {if(!is.na(subtitle)) ggplot2::labs(subtitle = subtitle)} +
    ggplot2::theme(legend.position = c(.75, .1),
                   legend.text = ggplot2::element_text(size = 15),
                   legend.direction = "horizontal", # legend.position = "none",
                   axis.title.x= ggplot2::element_blank(),
                   axis.text.x= ggplot2::element_blank(),
                   axis.ticks.x= ggplot2::element_blank())
  # Enrichment plot
  plt_e <- fgsea::plotEnrichment(
    pathway = geneset,
    stats = tibble::deframe(sig[,c("gene", rankcol)]) * -1 # flips mountain plt !!
  ) +
    ggplot2::ylab("Enrichment Score") +
    ggplot2::xlab("Rank") +
    ggplot2::theme_classic()

  grid <- cowplot::plot_grid(
    plt_wf, plt_e,
    align = "v", ncol = 1,
    rel_heights = c(2,1)
  )

  if(!is.null(savename)){
    cowplot::ggsave2(
      filename = savename,
      plot = grid,
      width = 9, height = 6
    )
  }

  return(grid)

}
