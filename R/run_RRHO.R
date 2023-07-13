
utils::globalVariables(c(
  "p1", "p2", "r1", "r2", "v1", "v2",
  ":=", "rank.1", "rank.2", "Unigene"
))

#' Run RRHO analysis
#'
#' Inference on the amount of agreement in two sorted lists using the Rank-Rank Hypergeometric Overlap test.
#' Outputs RRHO results, hypermatrix heatmap, rank rank scatter, metric scatter.
#'
#' Arrows in unicode are `\u2193` (down) and `\u2191` (up) if you want to be cool in your low/high descriptions B)
#'
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @encoding UTF-8
#'
#' @param sig1 string/dataframe; path or df for sig1, with cols `key` and `metric1`
#' @param sig2 string/dataframe; path or df for sig2, with cols `key` and `metric2`
#' @param sig1_name string; description of sig1
#' @param sig2_name string; description of sig2
#' @param sig1_low string; description of low sig1 values
#' @param sig1_high string; description of high sig1 values
#' @param sig2_low string; description of low sig2 values
#' @param sig2_high string; description of high sig2 values
#' @param key string; colname of corresponding values btwn both sigs
#' @param metric1 string; colname of sig1 metric to rank by
#' @param metric2 string; colname of sig2 metric to rank by
#' @param savename string; filepath to save results under (no extension)
#' @param webtool logical; T to output text formatted for Graeber RRHO webtool
#' @param steps integer; RRHO step size
#' @param BY logical; T for Benjamini-Yekutieli FDR corrected pvalues
#' @param hm_method string; make heatmap with ggplot geom_raster or lattice levelplot
#' @param palette string; RColorBrewer continuous palette name
#' @param waterfall logical; `TRUE` to include waterfall subplots (`hm_method = "ggplot"` only)
#' @param scatter logical; `TRUE` to output metric & rank scatterplot
#' @param plot_fmt string; file extension to save plot as
#'
#' @return RRHO results, as list of plots if `scatter == TRUE` or as ggplot RRHO heatmap object
#' @export

run_RRHO <- function(sig1, sig2, sig1_name, sig2_name,
                     sig1_low = NA, sig1_high = NA,
                     sig2_low = NA, sig2_high = NA,
                     key = "gene", metric1 = "sign_log_p", metric2 = metric1,
                     steps = NULL,
                     savename = NULL, webtool = TRUE, BY = FALSE,
                     hm_method = c("ggplot", "lattice"), palette = "Spectral",
                     waterfall = TRUE, scatter = TRUE, plot_fmt = "png"){
  hm_method <- match.arg(hm_method)

  # Individual signatures ----
  if(is.character(sig1)){sig1 <- Rubrary::rread(sig1)}
  m1 <- paste0(metric1, ".1")
  sig1 <- sig1 %>%
    select(!!sym(key), !!sym(metric1)) %>%
    filter(stats::complete.cases(.)) %>%
    distinct(!!sym(key), .keep_all = TRUE) %>%
    rename(!!m1 := !!sym(metric1))

  if(is.character(sig2)){sig2 <- Rubrary::rread(sig2)}
  m2 <- paste0(metric2, ".2")
  sig2 <- sig2 %>%
    select(!!sym(key), !!sym(metric2)) %>%
    filter(stats::complete.cases(.)) %>%
    distinct(!!sym(key), .keep_all = TRUE) %>%
    rename(!!m2 := !!sym(metric2))

  # Merge ----
  merged <- inner_join(x = sig1, y = sig2, by = key) %>%
    arrange(desc(!!sym(m2))) %>%
    mutate(rank.2 = 1:nrow(.)) %>%
    arrange(desc(!!sym(m1))) %>%
    mutate(rank.1 = 1:nrow(.)) %>%
    select(!!key, rank.1, rank.2, !!m1, !!m2)

  message(paste0("1) ", sig1_name, " genes: ", nrow(sig1)))
  message(paste0("2) ", sig2_name, " genes: ", nrow(sig2)))
  message(paste0("Intersect genes: ", nrow(merged)))

  if(is.null(steps)){
    defaultStepSize <- utils::getFromNamespace("defaultStepSize", "RRHO")
    steps <- defaultStepSize(list1 = sig1, list2 = sig2)
  }

  # Webtool ----
  if(webtool){
    df_web <- merged %>%
      rename(Unigene = !!sym(key),
             !!paste0("rank.", gsub(" ", "_", sig1_name)) := rank.1,
             !!paste0("rank.", gsub(" ", "_", sig2_name)) := rank.2,
             !!paste0(metric1,".", gsub(" ", "_", sig1_name)) := !!m1,
             !!paste0(metric2,".", gsub(" ", "_", sig2_name)) := !!m2) %>%
      mutate(Gene_Symbol = Unigene, .after = Unigene)

    Rubrary::rwrite(
      x = df_web,
      file = ifelse(is.null(savename), "RRHO_web.txt", paste0(savename, "_web.txt"))
    )

    # Output webtool inputs
    message("RRHO: https://systems.crump.ucla.edu/rankrank/rankranksimple.php")
    message(paste0("Dataset 1: ", sig1_name))
    if(!is.na(sig1_high) && !is.na(sig1_low)){
      message(paste0("** Class 1: ", sig1_low))
      message(paste0("** Class 2: ", sig1_high))
    }
    message(paste0("Dataset 2: ", sig2_name))
    if(!is.na(sig2_high) && !is.na(sig2_low)){
      message(paste0("** Class 1: ", sig2_low))
      message(paste0("** Class 2: ", sig2_high))
    }
    message(paste0("Step size: ", steps))
  }

  # Run RRHO ----
  Rubrary::use_pkg("RRHO")
  obj_RRHO <- RRHO::RRHO(
    list1 = sig1,
    list2 = sig2,
    labels = c(sig1_name, sig2_name),
    stepsize = steps,
    BY = BY,
    alternative = "enrichment"
  )

  # # Plot annotation with labels #n - ugly!!
  # sig1_low <- paste0(sig1_low, " (n = ", sum(sig1[,paste0(metric1,".1")] < 0),")")
  # sig1_high <- paste0(sig1_high, " (n = ", sum(sig1[,paste0(metric1,".1")] > 0),")")
  # sig2_low <- paste0(sig2_low, " (n = ", sum(sig2[,paste0(metric2,".2")] < 0),")")
  # sig2_high <- paste0(sig2_high, " (n = ", sum(sig2[,paste0(metric2,".2")] > 0),")")

  # Heatmap plot ----
  if(hm_method == "lattice"){ # lattice style heatmap
    requireNamespace("lattice", quietly = TRUE)
    plt <- lattice::levelplot(obj_RRHO$hypermat, # should be flipped?
                              xlab = '',
                              ylab = '',
                              col.regions=rev(grDevices::colorRampPalette(
                                RColorBrewer::brewer.pal(11, palette))(100)),
                              at = seq(0, 100, length = 100),
                              scales = list(tck = c(0,0),
                                            x=list(draw=FALSE),
                                            y=list(draw=FALSE)))

    if(!is.null(savename)){
      grDevices::png(filename = paste0(savename,"_heatmap.", plot_fmt),
          type = "cairo",
          units = "in",
          width = 5, height = 4, pointsize = 12, res = 96)
      print(plt)
      grDevices::dev.off()
    }
  } else { # ggplot style heatmap
    hypmtx_long <- obj_RRHO$hypermat %>%
      as.data.frame() %>%
      tibble::rownames_to_column("row") %>%
      tidyr::pivot_longer(!row, names_to = "col", values_to = "value") %>%
      mutate(col = sub("V", "", col)) %>%
      mutate_all(as.numeric)

    plt_hypmtx <- ggplot(hypmtx_long, aes(x = row, y = col, fill = value)) +
      geom_raster() +
      scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, palette)),
                           breaks = scales::pretty_breaks(5)) +
      theme_classic() +
      xlab(sig1_name) +
      ylab(sig2_name) +
      theme(line = element_blank(),
            axis.text = element_blank(),
            axis.title.x = element_text(size = rel(1.5), margin = margin(t = -10)),
            axis.title.y = element_text(size = rel(1.5), margin = margin(r = -10)),
            legend.title = element_blank(),
            legend.key.height = unit(0.125, "npc")) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      # scale_x_reverse(expand = c(0,0)) + # not sure if need to flip the raster...
      # scale_y_reverse(expand = c(0,0)) +
      # x-axis / sig1 annotations
      annotate("text", x = nrow(obj_RRHO$hypermat) + 6, y = -6, label = sig1_low, hjust = "right") +
      annotate("text", x = 0, y = -6, label = sig1_high, hjust = "left") +
      # y-axis / sig2 annotations
      annotate("text", x = -6, y = nrow(obj_RRHO$hypermat), label = sig2_low, hjust = "right", angle = 90) +
      annotate("text", x = -6, y = 0, label = sig2_high, hjust = "left", angle = 90) +
      coord_cartesian(clip = "off")

    ## Waterfall side panels ----
    if(waterfall){
      df_gg <- merged %>%
        `colnames<-`(c("id", "r1", "r2", "v1", "v2")) %>%
        mutate(p1 = ifelse(v1 >= 0, T, F),
               p2 = ifelse(v2 >= 0, T, F))
      # sig1 wf
      plt_wf1 <- ggplot(df_gg, aes(x = r1, y = v1, fill = p1)) +
        geom_col(show.legend = F) +
        scale_fill_manual(values = c("blue", "red")) +
        theme_classic() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 6)) +
        scale_x_continuous(breaks = scales::pretty_breaks(6), expand = c(0, 0), limits = c(NA, NA)) +
        scale_y_continuous(breaks = scales::pretty_breaks(6), expand = c(0, 0), limits = c(NA, NA))
      # sig2 wf
      plt_wf2 <- ggplot(df_gg, aes(x = r2, y = v2, fill = p2)) +
        geom_col(show.legend = F) +
        scale_fill_manual(values = c("blue", "red")) +
        theme_classic() +
        theme(axis.title = element_blank(),
              axis.text = element_text(size = 6)) +
        scale_x_continuous(breaks = scales::pretty_breaks(6), expand = c(0, 0), limits = c(NA, NA),
                           position = "top", guide = guide_axis(angle = 90)) +
        scale_y_reverse(breaks = scales::pretty_breaks(6), expand = c(0, 0), limits = c(NA, NA),
                        guide = guide_axis(angle = 90)) +
        coord_flip()

      grid_wf <- plt_wf2 + plt_hypmtx + plot_spacer() + plt_wf1 +
        plot_layout(widths = c(1, 4), heights = c(4, 1), guides = "keep")

      plt <- grid_wf
      wd <- 7; ht <- 6
    } else {
      plt <- plt_hypmtx
      wd <- 4; ht <- 3
    }
    if(!is.null(savename)){
      ggsave(
        filename = paste0(savename,"_heatmap.", plot_fmt),
        plot = plt,
        width = wd, height = ht
      )
    }
  }

  # Scatter plots ----
  if(scatter){
    ## Metric scatter ----
    metricsct <- Rubrary::plot_scatter(
      df = merged,
      xval = paste0(metric1, ".1"), xlabel = sig1_name,
      yval = paste0(metric2, ".2"), ylabel = sig2_name,
      title = paste0("Metric Scatter Plot"),
      pt_alpha = 0.5, reverse = F
    ) +
      scale_x_continuous(breaks = scales::pretty_breaks(6)) +
      scale_y_continuous(breaks = scales::pretty_breaks(6))

    ## Rank scatter ----
    ranksct <- Rubrary::plot_scatter(
      df = merged, rank = T, reverse = F,
      xval = paste0(metric1, ".1"), xlabel = sig1_name,
      yval = paste0(metric2, ".2"), ylabel = sig2_name,
      title = paste0("Rank Scatter Plot"),
      pt_alpha = 0.5, pt_size = 0.5
    ) +
      theme(line = element_blank(),
            axis.text = element_blank(),
            axis.title.x = element_text(size = rel(1.25), margin = margin(t = -(0.002 * nrow(merged)))),
            axis.title.y = element_text(size = rel(1.25), margin = margin(r = -(0.002 * nrow(merged))))) +
      geom_segment(x = 0, xend = nrow(merged), y = 0, yend = 0) +
      geom_segment(y = 0, yend = nrow(merged), x = 0, xend = 0) +
      # x-axis / sig1 annotations
      annotate("text", x = 0, y = -(0.025 * nrow(merged)), label = sig1_high,
               hjust = "left") +
      annotate("text", x = nrow(merged), y = -(0.025 * nrow(merged)), label = sig1_low,
               hjust = "right") +
      # y-axis / sig2 annotations
      annotate("text", x = -(0.025 * nrow(merged)), y = 0, label = sig2_high,
               hjust = "left", angle = 90) +
      annotate("text", x = -(0.025 * nrow(merged)), y = nrow(merged), label = sig2_low,
               hjust = "right", angle = 90) +
      coord_cartesian(clip = "off")

    if(waterfall){
      rnk_wf <- plt_wf2 + ranksct + plot_spacer() + plt_wf1 +
        plot_layout(widths = c(1, 4), heights = c(4, 1), guides = "keep")
      ranksct <- rnk_wf
    }

    if(!is.null(savename)){
      ggsave(
        filename = paste0(savename,"_metricsct.", plot_fmt),
        plot = metricsct,
        width = wd, height = wd # Square
      )

      ggsave(
        filename = paste0(savename,"_ranksct.", plot_fmt),
        plot = ranksct,
        width = wd, height = ht + 1
      )
    }
  }

  # Output ----
  if(scatter){
    RRHO <- list(
      heatmap = plt,
      metricsct = metricsct,
      ranksct = ranksct
    )
    class(RRHO) <- "RRHO"
  } else {
    RRHO <- plt
  }

  return(RRHO)
}
