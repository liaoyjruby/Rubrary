utils::globalVariables(c(
  ".data", "label_pos", "label", "rank_val"
))

#' Plot horizontal waterfall plot (with optional density)
#'
#' If `ggtext` package is installed, will do fancy version of p_enrichment. :)
#'
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#'
#' @param sig dataframe; names/genes in 1st column + `rankcol` columns
#' @param highlight vector; list of names to highlight in waterfall
#' @param rankcol string; colname of values to rank by
#' @param rankcol_name string; descriptive name of values to rank by
#' @param hightolow logical; T for x-axis by decreasing value
#' @param vertical logical; T for columns to be horizontal
#' @param label logical; T to label highlighted values
#' @param replab_size numeric; size of gene label
#' @param density logical; T for density plot as side panel
#' @param hllab string; description of highlighted values
#' @param otherlab string; description of non-highlighted values
#' @param pval logical; T to include KS enrichment pvalue
#' @param lab_high string; descriptor for high rankcol values
#' @param lab_low string; descriptor for low rankcol values
#' @param lab_size numeric; size of value annotation label
#' @param legendpos vector; value btwn 0-1 as legend coordinates (ggplot `legend.position`)
#' @param title string; plot title
#' @param colors vector; two colors, for highlight vs. other
#' @param width numeric; width of plot
#' @param height numeric; height of plot
#' @param savename string; filepath to save figure under
#'
#' @return Waterfall plot ggplot2 object
#' @export
#'
#' @examples
#' library(dplyr)
#' airway_deseq = Rubrary::airway_deseq_res %>% relocate(hgnc_symbol)
#' genes = Rubrary::GSEA_pathways$GOBP_REGULATION_OF_GLUCOSE_IMPORT[1:20]
#'
#' Rubrary::plot_waterfall(
#'   sig = airway_deseq,
#'   highlight = genes,
#'   rankcol = "sign_log_p",
#'   rankcol_name = "Sign log p value",
#'   lab_high = "\U2191 in treated\n\U2193 in untreated",
#'   lab_low = "\U2191 in untreated\n\U2193 in treated",
#'   title = "Airway Treated vs. Untreated DESeq - Glucose Import Regulation Genes",
#'   density = TRUE
#' )
#'
plot_waterfall <- function(
    sig, highlight, rankcol, rankcol_name = rankcol, hightolow = FALSE,
    vertical = FALSE, label = TRUE, replab_size = 3.5, density = FALSE,
    hllab = "Highlight", otherlab = "Other", pval = TRUE, colors = c("firebrick3", "gray"),
    lab_high = NULL, lab_low = NULL, lab_size = 6, legendpos = "none",
    title = NULL, savename = NULL, width = 10, height = 5){
  if(is.character(sig)){ sig <- Rubrary::rread(sig)}
  names(sig)[1] <- "name"
  if(hightolow){
    tmp = lab_high
    lab_high = lab_low
    lab_low = tmp
    sig <- sig %>% arrange(desc(!!sym(rankcol)))
  } else {
    sig <- sig %>% arrange(!!sym(rankcol))
  }
  sig <- sig %>%
    rename(rank_val = !!sym(rankcol)) %>%
    mutate(rank = 1:nrow(.),
           type = ifelse(name %in% highlight, hllab, otherlab),
           type = factor(type, levels = c(hllab, otherlab)))

  # Base plot ----
  wf <- ggplot(sig, aes(x = rank, y = rank_val)) +
    geom_bar(aes(fill = type), stat = "identity") +
    scale_fill_manual(values = colors) +
    xlab("Rank") +
    ylab(rankcol_name) +
    labs(title = title) +
    theme_classic()

  # Annotations calculations ----
  # Figure out dimensions and scale for label position
  layer <- layer_scales(wf)
  yrange <- layer$y$range$range # rankcol range
  ypos_low <- round((min(yrange) / 2), 2) # 4) * 3
  ypos_high <- round((max(yrange) / 4) * 3, 2)
  xrange <- layer$x$range$range # number of ranked items
  xpos_bot <- round(max(xrange) / 5, 2)
  xpos_mid <- round(max(xrange) / 2, 2)
  xpos_top <- round((max(xrange) / 5) * 4, 2)

  # Try to determine optimal label placement if rankcol both pos and neg
  # If more space above y = 0, then put label above y = 0
  # If more space below y = 0, put labels below y = 0
  # If about even, one above, one below?
  if (yrange[1] < 0) {
    yrng = yrange[2] - yrange[1]
    ydiff = abs(yrange[2]) - abs(yrange[1])
    if(abs(yrange[1]) < abs(yrange[2])){ # More space above
      ypos = ypos_high
    } else { # More space below
      ypos = ypos_low
    }
  }
  sig <- sig %>% mutate(label_pos = ifelse(rank_val < 0, 0, rank_val))
  nudge_lab <- max(sig$label_pos) / 10
  sig_hl <- sig %>% filter(type == hllab)

  ## Plot annotations ----
  if (vertical) {
    wf_lab <- wf +
      geom_segment( # Emphasize highlight genes by tracing over them
        data = sig_hl,
        aes(x = rank, xend = rank, y = 0, yend = rank_val), color = colors[1]) +
      coord_flip() +
      {if (label) ggrepel::geom_text_repel(
        data = sig_hl, aes(x = rank, y = label_pos, label = name),
        force = 2, hjust = 0, direction = "y", max.overlaps = Inf,
        size = replab_size, nudge_y = nudge_lab, segment.size = 0.1)} +
      {if (!is.null(lab_high)) geom_text(
        x = xpos_top, y = ypos_low, label = lab_high, size = lab_size)} +
      {if (!is.null(lab_low)) geom_text(
        x = xpos_bot, y = ypos_high, label = lab_low, size = lab_size)}
  } else {
    wf_lab <- wf +
      geom_segment(
        data = sig_hl, aes(x = rank, xend = rank, y = 0, yend = rank_val),
        color = colors[1]) +
      {if (label) ggrepel::geom_text_repel(
        data = sig_hl, aes(x = rank, y = label_pos, label = name),
        force = 2, angle = 90, hjust = 0, direction = "x", max.overlaps = Inf,
        size = replab_size, nudge_y = nudge_lab, segment.size = 0.1)} +
      {if (!is.null(lab_high)) geom_text(
        x = xpos_top, y = ypos, label = lab_high, size = lab_size)} +
      {if (!is.null(lab_low)) geom_text(
        x = xpos_bot, y = ypos, label = lab_low, size = lab_size)}
  }

  # KS p-value ----
  Rubrary::use_pkg("ggtext")
  if (!requireNamespace("ggtext", quietly = TRUE)){
    sbt <- "KS enrich. p-value = "
    plsbt <- element_text(size = 15)
  } else {
    sbt <- "p-val<sub>enrichment</sub> = "
    plsbt <- ggtext::element_markdown(size = 15)
  }
  wf_lab <- wf_lab +
    {if (pval) labs(subtitle = paste0( # "KS enrich. p-value = "
        sbt, signif(Rubrary::get_kspval(sig, "rank_val", "type", hllab), digits = 4)))} +
    theme(legend.title = element_blank(),
          legend.position = legendpos,
          plot.subtitle = plsbt)

  # Density side panel ----
  if (density) {
    dplt <- Rubrary::plot_density(
      df = sig[sig$type == hllab,],
      value = "rank",
      group = "type",
      xlab = "",
      colors = colors,
      pval = F) +
      scale_x_continuous(limits = c(1, max(sig$rank))) +
      theme(legend.position = "none",
            axis.text.y = element_blank())
    if (vertical) {
      dplt <- dplt + coord_flip()
      grid <- wf_lab + dplt +
        plot_layout(widths = c(3,1))
      width <- width + 3
    } else {
      grid <- wf_lab / dplt +
        plot_layout(heights = c(3, 1))
      height <- height + 3
    }
    plt <- grid
  } else {
    plt <- wf_lab
  }
  # Save ----
  if (!is.null(savename)) {
    ggsave(
      filename = savename,
      plot = plt,
      width = width, height = height
    )
  }
  return(plt)
}
