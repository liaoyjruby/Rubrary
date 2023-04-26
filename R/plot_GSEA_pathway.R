#' Plot GSEA pathway plots (waterfall + enrichment)
#'
#' @import ggplot2
#'
#' @param sig dataframe; "gene" + rankcol columns
#' @param geneset vector; list of genes
#' @param rankcol string; colname of values
#' @param rankcol_name string; descriptive name of values
#' @param hightolow logical; T for high values on left, low on right
#' @param label logical; T to label highlighted genes
#' @param legendpos vector; value btwn 0-1 as legend coordinates (ggplot legend.position option)
#' @param title string; plot title
#' @param subtitle string; plot subtitle; !! overwrites p-val_enrichment
#' @param savename string; filepath to save png under
#' @param highlab string; description of high values
#' @param lowlab string; description of low values
#' @param hllab string; description of highlighted values
#'
#' @return grid of gene set enrichment waterfall above enrichment plot
#' @export
#'
plot_GSEA_pathway <- function(sig, geneset, rankcol, rankcol_name = rankcol, hightolow = FALSE,
                              label = length(geneset) < 20, legendpos = "none",
                              highlab = NA, lowlab = NA, hllab = "Highlight",
                              title = "", subtitle = NULL, savename = NULL){
  Rubrary::use_pkg("fgsea")

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
    hightolow = hightolow
  ) +
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
    stats = tibble::deframe(sig[,c("gene", rankcol)])
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

#' `plot_GSEA_pathway` that works nicely with `lapply`
#'
#' @param path_name string; name of pathway
#' @param pthwys named list; key = geneset name, values = char vector of genes in geneset
#' @param sig_name string; name of signature
#' @param sig dataframe; signature
#' @param rankcol string; colname of values to rank by
#' @param rankcol_name string; descriptor of rankcol
#' @param hightolow logical; T for high values on left, low on right
#' @param label logical; T to label highlighted genes
#' @param hllab string; descriptor of highlighted genes
#' @param lowlab string; label for low rankcol values
#' @param highlab string; label for high rankcol values
#' @param savedir string; directory path for saving plot
#'
#' @return nothing; will save GSEA enrichment plots in savedir
#' @export
plot_GSEA_batch <- function(path_name, pthwys, sig_name, sig, rankcol, rankcol_name,
                            hllab = "Highlight", lowlab = "Low", highlab = "High",
                            savedir = "./", label = length(pthwys[[path_name]]) < 50,
                            hightolow = FALSE){
  Rubrary::plot_GSEA_pathway(
    sig = sig, legendpos = c(0.5, 0.8),
    rankcol = rankcol,
    rankcol_name = rankcol_name,
    geneset = pthwys[[path_name]],
    label = label,
    title = path_name,
    lowlab = lowlab,
    highlab = highlab,
    hllab = hllab,
    savename = paste0(savedir,"/",sig_name,"_", path_name, ".png")
  )
}
