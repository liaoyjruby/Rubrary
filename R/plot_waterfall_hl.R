#' Get number of TF MRA signature genes to highlight
#'
#' @param signature df; Has "pvalue" and "NES" columns
#' @param pos logical; "T" for positive NES
#' @param sig logical; "T" for p < 0.05
#'
#' @return integer with n number of genes to highlight
#' @export
get_n <- function(signature, pos = T, sig = F){
  if (sig) {
    signature <- signature[signature$pvalue < 0.05,]
  }
  if (pos){
    signature <- signature[signature$NES > 0,]
  } else {
    signature <- signature[signature$NES < 0,]
  }
  return(length(nrow(signature)))
}

#' Plot ranked waterfall with genes highlighted
#'
#' Revamp of Favour's original function
#'
#' @param sig df; Ranked signature with "rankcol" and "gene" columns
#' @param rankcol string; column to rank sig df by
#' @param sig_hl df; ranked signature with genes in column 1
#' @param topn integer; number of genes from sig_hl to highlight
#' @param rev logical; reverse to bottom n genes from sig_hl instead of top
#' @param toplabel string; name of highlight group
#' @param title string; plot title
#' @param subtitle string; plot subtitle
#' @param save logical; "T" to save
#' @param savename string; name to save plot under
#' @param label logical
#' @param ylab string
#' @param wf_pos string; label for positive side of waterfall
#' @param wf_neg string; label for negative side of waterfall
#'
#' @importFrom ggplot2 geom_text layer_scales
#' @importFrom ggpubr ggbarplot ggdensity rremove rotate
#'
#' @return Waterfall plot with n genes from sig_hl highlighted, ks p value
#' @export
plot_waterfall_hl <- function(sig, rankcol = "sign_log_p", sig_hl, topn = 200, rev = F, toplabel = "Top SCN",
                              label = F, ylab = "DESeq Signed log p-values",
                              wf_pos = "Post-Therapy", wf_neg = "Pre-Therapy",
                              title = "DE Genes", subtitle = NA,
                              save = F, savename = "DESeq_slogp_WF_hl.png") {

  # Rank DF by rankcol values
  sig <- sig[order(sig[,rankcol], decreasing = T),]
  sig$rank <- nrow(sig):1

  colors <- c("#120632", "#999999")
  if(!rev) {
    hlrange <- 1:topn
  } else {
    hlrange <- (nrow(sig_hl) - topn):nrow(sig_hl)
    # colors <- rev(colors)
  }
  sig$type <- ifelse(sig$gene %in% sig_hl[hlrange, 1], toplabel, "Other")
  sig$type <- factor(sig$type, levels = c(toplabel, "Other"))

  ks_pval = stats::ks.test(
    sig[sig$type == toplabel, "rank"],
    sig[!sig$type == toplabel, "rank"]
  )$p.value

  a = ggbarplot(
    data = sig,
    x = "rank",
    y = rankcol,
    xlab = "Gene Rank",
    ylab = ylab,
    title = title,
    rotate = T,
    palette = colors,
    fill = "type",
    color = "type",
    sort.val = "asc", # Ascending
    sort.by.groups = FALSE,
    legend.title = "Type",
    label = label, lab.pos = "out"
  )

  # Figure out dimensions and scale for label position
  layer <- layer_scales(a)
  yrange <- layer$y$range$range # slogp range
  ypos <- (yrange[which.max(abs(yrange))] / 2)
  xrange <- layer$x$range_c$range # number of genes
  xpos_bot <- max(xrange) / 4
  xpos_mid <- max(xrange) / 2
  xpos_top <- (max(xrange) / 4) * 3

  a <- a +
    rremove("y.ticks") +
    rremove("y.text") +
    theme(legend.position=c(0,1), legend.justification=c(0,1)) +
    geom_text(x = xpos_top, y = ypos, label = wf_pos) +
    geom_text(x = xpos_mid, y = ypos, label = paste0("KS enrich. p-value = ", signif(ks_pval, digits = 4))) +
    geom_text(x = xpos_bot, y = ypos, label = wf_neg) +
    {if(!is.na(subtitle)) labs(subtitle = subtitle)}

  x = ggdensity(
    data = sig,
    x = "rank",
    xlab = "",
    ylab = "Density",
    fill = "type",
    alpha = 0.85,
    palette = colors
  ) +
    rotate() +
    rremove("y.ticks") +
    rremove("y.text") +
    rremove("x.text") +
    rremove("legend")

  pg = cowplot::plot_grid(
    a,
    x,
    align = "h", axis = "bt",
    ncol = 2,
    rel_widths = c(2, 0.7), #increase 0.3 to 0.5 (or more) if density plot is too squished
    rel_heights = c(0.7, 2) #increase 0.3 to 0.5 (or more) if density plot is too squished
  )

  if(save){
    cowplot::ggsave2(
      filename = savename,
      plot = pg,
      width = 10, height = 8
    )
  }

  return(pg)
}
