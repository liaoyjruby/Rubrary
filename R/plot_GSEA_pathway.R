#' Plot GSEA pathway plots (waterfall + enrichment)
#'
#' @import ggplot2
#'
#' @param sig dataframe; "gene" + rankcol columns
#' @param geneset vector; list of genes
#' @param rankcol string; colname of values
#' @param rankcol_name string; descriptive name of values
#' @param label logical; T to label highlighted genes
#' @param legendpos vector; value btwn 0-1 as legend coordinates (ggplot legend.position option)
#' @param title string; title of plot
#' @param subtitle string; (name of signature?)
#' @param savename string; filepath to save png under
#' @param highlab string; description of high values
#' @param lowlab string; description of low values
#' @param hllab string; description of highlighted values
#' @param hightolow logical; left = high rankcol, right = low rankcol
#'
#' @return grid of gene set enrichment waterfall above enrichment plot
#' @export
#'
plot_GSEA_pathway <- function(sig, geneset, rankcol, rankcol_name = rankcol,
                              label = length(geneset) < 20, legendpos = "none",
                              hightolow = FALSE, highlab = NA, lowlab = NA,
                              hllab = "Highlight", title = "", subtitle = NULL, savename = NULL){
  path_genes <- geneset
  # Waterfall plot
  plt_wf <- Rubrary::plot_waterfall(
    sig = sig,
    label = label,
    highlight = path_genes,
    highlab = highlab,
    lowlab = lowlab,
    hllab = hllab,
    rankcol = rankcol,
    rankcol_name = rankcol_name,
    title = title,
  ) +
    {if(hightolow) scale_x_reverse()} +
    {if(!is.null(subtitle)) labs(subtitle = subtitle)} +
    theme(legend.position = legendpos,
                   legend.text = element_text(size = 15),
                   legend.direction = "horizontal", # legend.position = "none",
                   axis.title.x= element_blank(),
                   axis.text.x= element_blank(),
                   axis.ticks.x= element_blank())
  # Enrichment plot
  plt_e <- fgsea::plotEnrichment(
    pathway = geneset,
    stats = tibble::deframe(sig[,c("gene", rankcol)]) # flips mountain plt !!
  ) +
    ylab("Enrichment Score") +
    xlab("Rank") +
    theme_classic() +
    {if(!hightolow) scale_x_reverse()}

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
