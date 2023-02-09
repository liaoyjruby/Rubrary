utils::globalVariables(c(
  ".data", "label_pos", "label"
))

#' Plot waterfall plot (with optional density)
#'
#' @param sig dataframe
#' @param highlight vector
#' @param rankcol string
#' @param vert logical
#' @param density logical
#' @param ylab string
#' @param hllab string
#' @param otherlab string
#' @param pval logical
#' @param highlab string
#' @param lowlab string
#' @param title string
#' @param colors vector
#' @param width numeric
#' @param height numeric
#' @param savename string
#' @param label logical
#'
#' @return Waterfall plot
#' @importFrom ggplot2 ggplot aes geom_segment labs coord_flip layer_scales
#' @export
#'
plot_waterfall <- function(sig, highlight, rankcol, label = TRUE, vert = FALSE, density = FALSE,
                           ylab = rankcol, hllab = "Top SCN", otherlab = "Others",
                           pval = TRUE, highlab = NA, lowlab = NA,
                           title = NULL, colors = c("firebrick3", "gray"),
                           width = 12, height = 6,
                           savename = NA) {
  # Rank DF by rankcol values
  sig <- sig[order(sig[, rankcol], decreasing = F), ]
  sig$rank <- 1:nrow(sig)

  sig$type <- ifelse(sig[, 1] %in% highlight, hllab, otherlab)
  sig$type <- factor(sig$type, levels = c(hllab, otherlab))
  sig$label <- ifelse(sig[, 1] %in% highlight, sig[, 1], NA)
  sig_hl <- sig[sig[, 1] %in% highlight, ]

  # KS pval
  ks_pval <- stats::ks.test(
    sig[sig$type == hllab, "rank"],
    sig[!sig$type == hllab, "rank"]
  )$p.value

  # Base barplot
  wf <- ggpubr::ggbarplot(
    data = sig,
    x = "rank",
    y = rankcol,
    xlab = "Rank",
    ylab = ylab,
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
  ypos_left <- (min(yrange) / 4) * 3
  ypos_right <- (max(yrange) / 4) * 3
  xrange <- layer$x$range$range # number of ranked items
  xpos_bot <- max(xrange) / 4
  xpos_mid <- max(xrange) / 2
  xpos_top <- (max(xrange) / 4) * 3

  sig$label_pos <- ifelse(sig[, rankcol] < 0, 0, sig[, rankcol])
  nudge_lab <- max(sig$label_pos) / 10

  if (vert) {
    wf_lab <- wf +
      geom_segment(data = sig_hl, aes(x = rank, xend = rank, y = 0, yend = .data[[rankcol]], ), color = colors[1]) +
      coord_flip() +
      {if (!is.na(highlab)) geom_text(x = xpos_top, y = ypos_left, label = highlab)} +
      {if (!is.na(lowlab)) geom_text(x = xpos_bot, y = ypos_right, label = lowlab)} +
      {if (label) ggrepel::geom_text_repel(
        data = sig,
        aes(x = rank, y = label_pos, label = label),
        force = 2, hjust = 0, direction = "y",
        size = 4, nudge_y = nudge_lab, segment.size = 0.1
      )}
  } else {
    wf_lab <- wf +
      {if (label) ggrepel::geom_text_repel(
        data = sig,
        aes(x = rank, y = label_pos, label = label),
        force = 2, angle = 90, hjust = 0, direction = "x",
        size = 4, nudge_y = nudge_lab, segment.size = 0.1
      )} +
      geom_segment(
        data = sig_hl,
        aes(x = rank, xend = rank, y = 0, yend = .data[[rankcol]], ),
        color = colors[1]
      ) +
      {if (!is.na(highlab)) geom_text(x = xpos_top, y = yrange[2], label = highlab)} +
      {if (!is.na(lowlab)) geom_text(x = xpos_bot, y = yrange[2], label = lowlab)}
  }

  # Add pvalue, remove legend title
  wf_lab <- wf_lab +
    {if (pval) labs(subtitle = paste0("KS enrich. p-value = ", signif(ks_pval, digits = 4)))} +
    theme(legend.title = element_blank())

  # Save
  if (!is.na(savename)) {
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
