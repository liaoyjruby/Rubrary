#' Plot GSEA pathway plots (waterfall + enrichment)
#'
#' If subtitle is not specified, Kolmogorov-Smirnov enrichment p-value is calculated to check if the genes in the pathway separate from the genes not in the chosen pathway based on ranking by the desired `rankcol`.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param sig dataframe; genecol + rankcol columns
#' @param geneset vector; list of genes
#' @param genecol string; colname of gene names in `sig`
#' @param rankcol string; colname of values in `sig`
#' @param rankcol_name string; descriptive name of values
#' @param hightolow logical; T for high values on left, low on right
#' @param label logical; T to label highlighted genes
#' @param legendpos vector; value btwn 0-1 as legend coordinates (ggplot2 legend.position)
#' @param title string; plot title
#' @param subtitle string; plot subtitle; !! overwrites p-val_enrichment
#' @param savename string; filepath to save png under
#' @param lab_high string; description of high values
#' @param lab_low string; description of low values
#' @param hllab string; description of highlighted values
#'
#' @return grid of gene set enrichment waterfall above enrichment plot
#' @export
#'
plot_GSEA_pathway <- function(sig, geneset, genecol = "gene", rankcol, rankcol_name = rankcol, hightolow = FALSE,
                              label = length(geneset) < 20, legendpos = "none",
                              lab_high = NULL, lab_low = NULL, hllab = "Highlight",
                              title = NULL, subtitle = NULL, savename = NULL){
  Rubrary::use_pkg("fgsea")

  sig <- sig %>%
    rename(gene = any_of(genecol)) %>%
    select(gene, everything())

  path_genes <- geneset
  # Waterfall plot
  plt_wf <- Rubrary::plot_waterfall(
    sig = sig,
    label = label,
    highlight = path_genes,
    lab_high = lab_high,
    lab_low = lab_low,
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
#' @param sig dataframe; signature
#' @param genecol string; colname of gene names in `sig`
#' @param rankcol string; colname of values to rank by
#' @param rankcol_name string; descriptor of rankcol
#' @param hllab string; descriptor of highlighted genes
#' @param hightolow logical; T for high values on left, low on right
#' @param lab_low string; label for low rankcol values
#' @param lab_high string; label for high rankcol values
#' @param legendpos vector; value btwn 0-1 as legend coordinates (ggplot2 legend.position)
#' @param label logical; T to label points in plot
#' @param sig_name string; name of signature
#' @param savedir string; directory path for saving plot
#'
#' @return Gene set enrichment plot as ggplot object
#' @export
plot_GSEA_batch <- function(path_name, pthwys, sig, genecol = "gene", rankcol, rankcol_name = rankcol,
                            hllab = "Pathway genes", hightolow = FALSE,
                            lab_low = NULL, lab_high = NULL, legendpos = c(0.5, 0.2),
                            label = length(pthwys[[path_name]]) < 50,
                            sig_name = "", savedir = NULL){
  if(!is.null(savedir)){
    sig_name <- if(sig_name != "") paste0(sig_name,"_")
    savename = paste0(savedir,"/",sig_name, path_name, ".png")
  } else {
    savename = NULL
  }

  Rubrary::plot_GSEA_pathway(
    sig = sig, legendpos = legendpos,
    genecol = genecol,
    rankcol = rankcol,
    rankcol_name = rankcol_name,
    geneset = pthwys[[path_name]],
    label = label,
    title = path_name,
    lab_low = lab_low,
    lab_high = lab_high,
    hllab = hllab,
    savename = savename
  )
}
