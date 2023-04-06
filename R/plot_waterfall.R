utils::globalVariables(c(
  ".data", "label_pos", "label"
))

#' Plot waterfall plot (with optional density)
#'
#' @param sig dataframe; name (1st column) + rankcol columns
#' @param highlight vector; list of names to highlight in waterfall
#' @param rankcol string; colname of values
#' @param rankcol_name string; descriptive name of values
#' @param label logical; T to label highlighted values
#' @param vert logical; T for columns to be horizontal (not optimized)
#' @param density logical; T to output aligned density plot w/ wf
#' @param hllab string; description of highlighted values
#' @param otherlab string; description of non-highlighted values
#' @param pval logical; T to include KS enrichment pvalue
#' @param highlab string; descriptor for high rankcol values
#' @param lowlab string; descriptor for low rankcol values
#' @param legendpos vector; value btwn 0-1 as legend coordinates (ggplot legend.position option)
#' @param title string; plot title
#' @param colors vector; two colors, for highlight vs. other
#' @param width numeric; width of plot
#' @param height numeric; height of plot
#' @param savename string; filepath to save figure under
#'
#' @return Waterfall plot
#' @importFrom ggplot2 ggplot aes geom_segment labs coord_flip layer_scales element_text
#' @export
#'
plot_waterfall <- function(sig, highlight, rankcol, rankcol_name = rankcol, label = TRUE,
                           vert = FALSE, density = FALSE, hllab = "Highlight", otherlab = "Others",
                           pval = TRUE, highlab = NA, lowlab = NA, legendpos = "none",
                           title = NULL, colors = c("firebrick3", "gray"),
                           width = 10, height = 5,
                           savename = NULL) {
  # Rank DF by rankcol values
  sig <- sig[order(sig[, rankcol], decreasing = F), ]
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
  ypos <- (yrange[which.max(abs(yrange))] / 2)
  ypos_left <- (min(yrange) / 2) # 4) * 3
  ypos_right <- (max(yrange) / 4) * 3
  xrange <- layer$x$range$range # number of ranked items
  xpos_bot <- max(xrange) / 5
  xpos_mid <- max(xrange) / 2
  xpos_top <- (max(xrange) / 5) * 4
  lab_size <- 6
  replab_size <- 3.5

  sig$label_pos <- ifelse(sig[, rankcol] < 0, 0, sig[, rankcol])
  nudge_lab <- max(sig$label_pos) / 10

  # Place high/low text annotations
  if (vert) {
    wf_lab <- wf +
      geom_segment(data = sig_hl, aes(x = rank, xend = rank, y = 0, yend = .data[[rankcol]], ), color = colors[1]) +
      coord_flip() +
      {if (!is.na(highlab)) geom_text(x = xpos_top, y = ypos_left, label = highlab, size = lab_size)} +
      {if (!is.na(lowlab)) geom_text(x = xpos_bot, y = ypos_right, label = lowlab, size = lab_size)} +
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
      {if (!is.na(highlab)) geom_text(x = xpos_top, y = ypos_left, label = highlab, size=lab_size)} +
      {if (!is.na(lowlab)) geom_text(x = xpos_bot, y = ypos_left, label = lowlab, size=lab_size)}
  }

  # Add pvalue, remove legend title
  wf_lab <- wf_lab +
    {if (pval) labs(subtitle = paste0("p-val<sub>enrichment</sub> = ",
                                      signif(Rubrary::get_kspval(sig, rankcol, "type", hllab), digits = 4)))} + # "KS enrich. p-value = "
    theme(legend.title = element_blank(),
          legend.position = legendpos,
          plot.subtitle = ggtext::element_markdown(size = 15))

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
