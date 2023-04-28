utils::globalVariables(c(
  ".data", "label_pos", "label"
))

#' Plot horizontal waterfall plot (with optional density)
#'
#' If `ggtext` package is installed, will do fancy version of p_enrichment. :)
#'
#' @import ggplot2
#'
#' @param sig dataframe; name (1st column) + rankcol columns
#' @param highlight vector; list of names to highlight in waterfall
#' @param rankcol string; colname of values
#' @param rankcol_name string; descriptive name of values
#' @param hightolow logical; T for x-axis by decreasing value
#' @param vert (NOT TESTED) logical; T for columns to be horizontal
#' @param label logical; T to label highlighted values
#' @param replab_size numeric; size of gene label
#' @param density logical; T to output aligned density plot w/ wf
#' @param hllab string; description of highlighted values
#' @param otherlab string; description of non-highlighted values
#' @param pval logical; T to include KS enrichment pvalue
#' @param lab_high string; descriptor for high rankcol values
#' @param lab_low string; descriptor for low rankcol values
#' @param lab_size numeric; size of value annotation label
#' @param legendpos vector; value btwn 0-1 as legend coordinates (ggplot legend.position option)
#' @param title string; plot title
#' @param colors vector; two colors, for highlight vs. other
#' @param width numeric; width of plot
#' @param height numeric; height of plot
#' @param savename string; filepath to save figure under
#'
#' @return Waterfall plot ggplot2 object
#' @export
#'
plot_waterfall <- function(sig, highlight, rankcol, rankcol_name = rankcol, hightolow = FALSE,
                           vert = FALSE, label = TRUE, replab_size = 3.5, density = FALSE,
                           hllab = "Highlight", otherlab = "Others", pval = TRUE,
                           lab_high = NULL, lab_low = NULL, lab_size = 6, legendpos = "none",
                           title = NULL, colors = c("firebrick3", "gray"),
                           width = 10, height = 5,
                           savename = NULL) {
  if(hightolow){
    tmp = lab_high
    lab_high = lab_low
    lab_low = tmp
  }

  # Rank DF by rankcol values
  sig <- sig[order(sig[, rankcol], decreasing = hightolow), ]
  sig$rank <- 1:nrow(sig)

  sig$type <- ifelse(sig[, 1] %in% highlight, hllab, otherlab)
  sig$type <- factor(sig$type, levels = c(hllab, otherlab))
  sig$label <- ifelse(sig[, 1] %in% highlight, sig[, 1], NA)
  sig_hl <- sig[sig[, 1] %in% highlight, ]

  # Base barplot
  wf <- ggpubr::ggbarplot(
    data = sig,
    x = "rank",
    y = rankcol,
    xlab = "Rank",
    ylab = rankcol_name,
    title = title,
    palette = colors,
    fill = "type",
    color = "type"
  )

  # Additional annotations
  # Figure out dimensions and scale for label position
  layer <- layer_scales(wf)
  yrange <- layer$y$range$range # rankcol range
  ypos_low <- (min(yrange) / 2) # 4) * 3
  ypos_high <- (max(yrange) / 4) * 3
  xrange <- layer$x$range$range # number of ranked items
  xpos_bot <- max(xrange) / 5
  xpos_mid <- max(xrange) / 2
  xpos_top <- (max(xrange) / 5) * 4

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


  sig$label_pos <- ifelse(sig[, rankcol] < 0, 0, sig[, rankcol])
  nudge_lab <- max(sig$label_pos) / 10

  # Place high/low text annotations
  if (vert) {
    wf_lab <- wf +
      geom_segment(data = sig_hl, aes(x = rank, xend = rank, y = 0, yend = .data[[rankcol]], ), color = colors[1]) +
      coord_flip() +
      {if (!is.null(lab_high)) geom_text(x = xpos_top, y = ypos_low, label = lab_high, size = lab_size)} +
      {if (!is.null(lab_low)) geom_text(x = xpos_bot, y = ypos_high, label = lab_low, size = lab_size)} +
      {if (label) ggrepel::geom_text_repel(
        data = sig,
        aes(x = rank, y = label_pos, label = label),
        force = 2, hjust = 0, direction = "y",
        size = replab_size, nudge_y = nudge_lab, segment.size = 0.1
      )}
  } else {
    wf_lab <- wf +
      {if (label) ggrepel::geom_text_repel(
        data = sig,
        aes(x = rank, y = label_pos, label = label),
        force = 2, angle = 90, hjust = 0, direction = "x",
        size = replab_size, nudge_y = nudge_lab, segment.size = 0.1
      )} +
      geom_segment(
        data = sig_hl,
        aes(x = rank, xend = rank, y = 0, yend = .data[[rankcol]], ),
        color = colors[1]
      ) +
      {if (!is.null(lab_high)) geom_text(x = xpos_top, y = ypos, label = lab_high, size=lab_size)} +
      {if (!is.null(lab_low)) geom_text(x = xpos_bot, y = ypos, label = lab_low, size=lab_size)}
  }

  if (!requireNamespace("ggtext", quietly = TRUE)){
    sbt <- "KS enrich. p-value = "
    plsbt <- element_text(size = 15)
  } else {
    sbt <- "p-val<sub>enrichment</sub> = "
    plsbt <- ggtext::element_markdown(size = 15)
  }

  # Add pvalue, remove legend title
  wf_lab <- wf_lab +
    {if (pval) labs(subtitle = paste0(sbt,
                                      signif(Rubrary::get_kspval(sig, rankcol, "type", hllab), digits = 4)))} + # "KS enrich. p-value = "
    theme(legend.title = element_blank(),
          legend.position = legendpos,
          plot.subtitle = plsbt)

  # Save
  if (!is.null(savename)) {
    ggsave(
      filename = savename,
      plot = wf_lab,
      width = width, height = height
    )
  }

  if (density) {
    dplt <- plot_density(
      df = sig,
      value = "rank",
      group = "type",
      xlab = "",
      colors = colors,
      pval = F
      ) +
      theme(legend.position = "none", axis.text.y = element_blank())

    if (vert) {
      dplt <- dplt + coord_flip()
      grid <- cowplot::plot_grid(
        wf_lab, dplt,
        align = 'h', nrow = 1,
        rel_widths = c(3,1)
      )
      width = width + 3

    } else {
      grid <- cowplot::plot_grid(
        wf_lab, dplt,
        align = 'v', ncol = 1,
        rel_heights = c(3,1)
      )
      height = height + 3

    }

    cowplot::ggsave2(
      filename = paste0(tools::file_path_sans_ext(savename), "_density_grid.png"),
      plot = grid,
      width = width,
      height = height
    )
  }

  return(wf_lab)
}
