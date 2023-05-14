utils::globalVariables(c(
  "NES", "name", "ord", "pathway", "pos", "ticks", "zero"
))

#' Format GSEA pathway name
#'
#' Strips prepended database abbreviation, changes underscores to spaces, and converts to title case.
#'
#' @param pw string; MSigDB pathway name
#' @param ignore char vector; terms to leave as all uppercase (ex. DNA, RNA) instead of title case
#' @param split logical; (roughly) split long names into multiple line?
#' @param split_nchar integer; if `split = TRUE`, max number of characters per line
#'
#' @return Nicer looking title-case'd pathway name
#' @export
#'
#' @examples
#' format_GSEA_name("DATABASE_DNA_REPAIR")
#' format_GSEA_name("DATABASE_DNA_REPAIR_THIS_IS_A_REALLY_LONG_PATHWAY", split = TRUE)
format_GSEA_name <- function(pw, ignore = c("DNA", "RNA", "NADH", " IV "),
                             split = FALSE, split_nchar = 40){
  Rubrary::use_pkg("stringr")
  rep <- stats::setNames(toupper(ignore), tolower(ignore))
  name <- sub("(.*?)[_(.*?)]", "", pw) %>% # Remove GSEA source prefix
    gsub("_", " ", .) %>%
    tolower() %>%
    stringr::str_replace_all(rep) %>% # Like a vectorized gsub
    tools::toTitleCase()
  if(split){ name <- lapply(name, Rubrary::split_line, chars = split_nchar) }
  return(name)
}

#' Plot GSEA barplot
#'
#' Can plot NES of pathways from one GSEA result, or compare NES for chosen pathways between two GSEA results.
#'
#' @import ggplot2
#'
#' @param gsea_res dataframe; GSEA results w/ `pathway` and `NES` columns
#' @param gsea_name string; description of GSEA results
#' @param gsea_pws char vector; pathways in GSEA results to plot
#' @param n_pws integer; if no pathways provided, top/bottom n pathways to plot
#' @param pw_format logical; clean up pathway names?
#' @param pw_split logical; split pathway names into multiple lines?
#' @param pw_size numeric; pathway name text size
#' @param gsea2_res dataframe; 2nd GSEA results w/ `pathway` and `NES` columns
#' @param gsea2_name string; description of 2nd GSEA results
#' @param order2 logical; order pathways by NES of 2nd GSEA results?
#' @param ptrn2 string; `ggpattern` pattern arg: 'stripe', 'crosshatch', 'point', 'circle', 'none'
#' @param NES_cutoff numeric; value to draw NES cutoff line at
#' @param cols vector; `cols[1]` for positive bar color, `cols[2]` for negative bar color
#' @param title string; plot title
#' @param savename string; file path to save plot under
#' @param width numeric; plot width
#' @param height numeric; plot height
#'
#' @return GSEA NES barplot as `ggplot` object
#' @export
#'
plot_GSEA_barplot <- function(gsea_res, gsea_name = "GSEA", gsea_pws = NULL, n_pws = 5,
                              pw_format = FALSE, pw_split = FALSE, pw_size = 5,
                              gsea2_res = NULL, gsea2_name = NULL, order2 = FALSE, ptrn2 = "stripe",
                              NES_cutoff = 1, cols = c("firebrick", "darkblue"), title = NULL,
                              savename = NULL, width = 8, height = NULL){
  #### CLEAN
  if(is.null(gsea_pws)){
    gsea_pws <- c(gsea_res[order(gsea_res$NES, decreasing = T),]$pathway[1:n_pws],
                  gsea_res[order(gsea_res$NES, decreasing = F),]$pathway[1:n_pws])
    n_pws <- 2*n_pws
  }

  GSEA <- gsea_res %>%
    filter(pathway %in% gsea_pws) %>%
    arrange(desc(NES)) %>%
    mutate(name = gsea_name)

  pw_order <- factor(GSEA$pathway, levels = GSEA$pathway)
  n_pos <- nrow(GSEA[GSEA$NES > 0,])
  n_pws <- length(GSEA$pathway)
  n_res <- 1

  if(!is.null(gsea2_res)){ # 2nd GSEA results
    Rubrary::use_pkg("ggpattern")
    gsea2_res <- gsea2_res %>%
      filter(pathway %in% gsea_pws) %>%
      arrange(desc(NES)) %>%
      mutate(name = gsea2_name)
    if(order2){ # Order by 2nd GSEA results NES instead
      pw_order <- factor(gsea2_res$pathway, levels = gsea2_res$pathway)
      n_pos <- nrow(gsea2_res[gsea2_res$NES > 0,])
    }
    GSEA <- rbind(GSEA, gsea2_res)
    n_res = 2
  }

  GSEA <- GSEA %>%
    mutate(pos = factor(ifelse(NES > 0, "Pos", "Neg"), levels = c("Pos", "Neg")),
           name = factor(name, levels = c(gsea2_name, gsea_name)),
           pathway = factor(pathway, levels = pw_order)) %>%
    arrange(pathway) %>%
    # "Pseudo" discrete scale - axis separates pos/neg NES pws so leave a gap
    mutate(ord = c(rep(1:n_pos, each = n_res),
                   rep((n_pos + 2):(length(unique(pathway))+1), each = n_res)))
  if(pw_format){
    GSEA$pw_name <- format_GSEA_name(GSEA$pathway, split = pw_split)
  } else {
    GSEA$pw_name <- GSEA$pathway
  }
  pw_labs <- c(GSEA$pw_name[1:n_pos], "", GSEA$pw_name[(n_pos+1):length(pw_order)])

  #### PLOT
  pseud_0 = n_pos + 0.75 # Midpoint in pseudofactor
  pseud_ord = c(1:n_pos, pseud_0, (n_pos+2):(length(pw_order)+1)) # Pseudo-order /w midpoint
  # Tick frame - infer from plot?
  # https://stackoverflow.com/questions/17753101/center-x-and-y-axis-with-ggplot2
  tick_frame <- data.frame(ticks = seq(-2, 2, length.out = 5), zero = pseud_0)

  plt <- ggplot(GSEA, aes(x = ord, y = NES)) +
    xlab("NES") +
    #### AXES
    geom_vline(xintercept = pseud_0) + # NES / x-axis
    geom_segment(aes(x = 0.25, xend = pseud_0, # Pathway / y-axis, pos side
                     y = 0, yend = 0)) +
    geom_segment(aes(y = 0, yend = 0, # Pathway / y-axis, neg side
                     x = pseud_0 + 0.5, xend = max(pseud_ord) + 0.75)) +
    scale_x_reverse(breaks = pseud_ord,
                    labels = pw_labs) +
    #### TICKS
    ## NES
    geom_segment(data = tick_frame,
                 aes(x = zero, xend = zero + .1,
                     y = ticks, yend = ticks)) +
    ## Pathway
    # NES positive
    geom_segment(data = GSEA[GSEA$pos == "Pos",],
                 aes(x = ord, xend = ord, y = 0, yend = -0.025)) +
    # NES negative
    geom_segment(data = GSEA[GSEA$pos == "Neg",],
                 aes(x = ord, xend = ord, y = 0, yend = 0.025)) +
    # NES tick values
    geom_label(data = tick_frame, aes(x = zero, y = ticks, label = ticks),
               vjust = 1.75, size = 3.5, label.size = 0., label.padding = unit(0.1, "lines")) +
    # Pathway names next to ticks - make conditional based on if paired pws have same sign?
    # NES positive
    geom_text(data = GSEA[GSEA$pos == "Pos",], aes(x = ord, y = 0),
              label = GSEA[GSEA$pos == "Pos",]$pw_name, size = pw_size,
              hjust = 1, nudge_y = -0.075) +
    # NES negative
    geom_text(data = GSEA[GSEA$pos == "Neg",], aes(x = ord, y = 0),
              label = GSEA[GSEA$pos == "Neg",]$pw_name, size = pw_size,
              hjust = 0, nudge_y = 0.075)

  #### BARPLOT
  if(!is.null(gsea2_res)){ # 2nd GSEA
    plt <- plt +
      ggpattern::geom_bar_pattern(aes(x = ord, y = NES, group = name, pattern = name, fill = pos),
                                  stat = "identity", position = position_dodge(), color = "black",
                                  pattern_fill = "white", pattern_color = "white",
                                  pattern_spacing = 0.01) +
      ggpattern::scale_pattern_manual(values = c(ptrn2, "none")) +
      ggpattern::scale_pattern_color_manual(values = cols[1])
  } else {
    plt <- plt +
      geom_bar(aes(fill = pos),
               stat = "identity", position = position_dodge(), color = "black")
  }

  plt <- plt +
    scale_fill_manual(values = cols)
  #### NES = +-1 guidelines - make conditional segment/hline if paired not corresponding?
  plt <- plt +
    geom_segment(aes(x = 0.5, xend = pseud_0 - 0.25,
                     y = NES_cutoff, yend = NES_cutoff), linetype = "dashed") +
    geom_segment(aes(x = pseud_0 + 0.75, xend = max(pseud_ord) + 0.5,
                     y = -1 * NES_cutoff, yend = -1 * NES_cutoff), linetype = "dashed") +
    # geom_hline(yintercept = -1, linetype = "dashed") +

    #### MISC
    guides(fill = "none") +
    labs(
      title = title,
    ) +
    theme_classic() + coord_flip() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 15),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.525, size = 15), # "NES" annotation
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())

  if(!is.null(savename)){
    if(is.null(height)){ height = n_pws }
    ggsave(
      filename = savename,
      plot = plt,
      width = width, height = height
    )
  }
  return(plt)
}

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
plot_GSEA_pathway_batch <- function(path_name, pthwys, sig, genecol = "gene", rankcol, rankcol_name = rankcol,
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
