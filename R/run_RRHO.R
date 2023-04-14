
utils::globalVariables(c(
  "p1", "p2", "r1", "r2", "v1", "v2"
))

#' Run RRHO analysis
#'
#' Inference on the amount of agreement in two sorted lists using the Rank-Rank Hypergeometric Overlap test.
#' Outputs RRHO results, hypermatrix heatmap, rank rank scatter, metric scatter
#'
#' @import ggplot2
#' @import patchwork
#'
#' @param sig1 string/dataframe; path or df for sig1, with cols "key" and "metric1"
#' @param sig2 string/dataframe; path or df for sig2, with cols "key" and "metric2"
#' @param sig1_name string; description of sig1
#' @param sig2_name string; description of sig2
#' @param sig1_low string; description of low sig1 values
#' @param sig1_high string; description of high sig1 values
#' @param sig2_low string; description of low sig2 values
#' @param sig2_high string; description of high sig2 values
#' @param key string; colname of corresponding values btwn both sigs
#' @param metric1 string; colname of sig1 metric to rank by
#' @param metric2 string; colname of sig2 metric to rank by
#' @param savename string; filepath to save results under (no ext.)
#' @param webtool logical; T to output text formatted for Graeber RRHO webtool
#' @param steps integer; RRHO step size
#' @param BY logical; T for Benjamini-Yekutieli FDR corrected pvalues
#' @param hm_method string; make heatmap with ggplot geom_raster or lattice levelplot
#' @param palette string; RColorBrewer continuous palette name
#' @param waterfall logical; T to include waterfall subplots (plot_method = "ggplot" only)
#' @param scatter logical; T to output metric & rank scatterplot
#'
#' @return RRHO results
#' @export

run_RRHO <- function(sig1, sig2, sig1_name, sig2_name,
                     sig1_low = NA, sig1_high = NA,
                     sig2_low = NA, sig2_high = NA,
                     key = "gene", metric1 = "sign_log_p", metric2 = metric1,
                     steps = NULL,
                     savename = NULL, webtool = TRUE, BY = FALSE,
                     hm_method = c("ggplot", "lattice"), palette = "Spectral",
                     waterfall = TRUE, scatter = TRUE){
  Rubrary::use_pkg("RRHO")

  if(is.null(steps)){
    defaultStepSize <- utils::getFromNamespace("defaultStepSize", "RRHO")
  }

  hm_method <- match.arg(hm_method)

  # Load & clean individual sigs
  if(is.character(sig1)){sig1 <- utils::read.delim(sig1)}
  sig1 <- sig1[,c(key, metric1)]
  sig1 <- sig1[stats::complete.cases(sig1),]
  sig1 <- sig1[!duplicated(sig1[,key]),]
  colnames(sig1) <- c(key, paste0(metric1,".1"))

  if(is.character(sig2)){sig2 <- utils::read.delim(sig2)}
  sig2 <- sig2[,c(key, metric2)]
  sig2 <- sig2[stats::complete.cases(sig2),]
  sig2 <- sig2[!duplicated(sig2[,key]),]
  colnames(sig2) <- c(key, paste0(metric2,".2"))

  # Join and order
  merged <- dplyr::inner_join(x = sig1, y = sig2, by = key)
  merged <- merged[order(merged[,paste0(metric1,".1")], decreasing = TRUE),]
  merged$rank.1 <- 1:nrow(merged)
  merged <- merged[order(merged[,paste0(metric2,".2")], decreasing = TRUE),]
  merged$rank.2 <- 1:nrow(merged)
  merged <- merged[c(key, "rank.1", "rank.2", paste0(metric1,".1"), paste0(metric1,".2"))]
  merged <- merged[order(merged$rank.1),]

  message(paste0("1) ", sig1_name, " genes: ", nrow(sig1)))
  message(paste0("2) ", sig2_name, " genes: ", nrow(sig2)))
  message(paste0("Intersect genes: ", nrow(merged)))

  if(webtool){
    df_web <- data.frame(
      Unigene = merged[,1],
      Gene_Symbol = merged[,1]
    )
    df_web[, paste0("rank.", gsub(" ", "_", sig1_name))] <- merged$rank.1
    df_web[, paste0("rank.", gsub(" ", "_", sig2_name))] <- merged$rank.2
    df_web[, paste0(metric1,".", gsub(" ", "_", sig1_name))] <- merged[,paste0(metric1,".1")]
    df_web[, paste0(metric2,".", gsub(" ", "_", sig2_name))] <- merged[,paste0(metric2,".2")]

    write.table(
      x = df_web,
      file = ifelse(is.null(savename), "RRHO_web.txt", paste0(savename, "_web.txt")),
      sep = "\t", quote = F, row.names = F
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

  # Run RRHO in R
  obj_RRHO <- RRHO::RRHO(
    list1 = sig1,
    list2 = sig2,
    labels = c(sig1_name, sig2_name),
    stepsize = steps,
    BY = BY,
    alternative = "enrichment"
  )

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
      grDevices::png(filename = paste0(savename,"_heatmap.png"),
          type = "cairo",
          units = "in",
          width = 5, height = 4, pointsize = 12, res = 96)
      print(plt)
      grDevices::dev.off()
    }
  } else { # ggplot style heatmap
    # Format data to long
    hypmtx <- tibble::rownames_to_column(as.data.frame(obj_RRHO$hypermat), var = "row")
    hypmtx_long <- reshape2::melt(
      data = hypmtx,
      id.vars = "row",
      variable.name = "col"
    )
    hypmtx_long$col <- sub("V", "", hypmtx_long$col)
    hypmtx_long <- dplyr::mutate_all(hypmtx_long, function(x) as.numeric(as.character(x)))

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
      scale_x_reverse(expand = c(0,0)) +
      scale_y_reverse(expand = c(0,0)) +
      # x-axis / sig1 annotations
      annotate("text", x = 0, y = nrow(hypmtx) + 6, label = sig1_low, hjust = "right") +
      annotate("text", x = nrow(hypmtx), y = nrow(hypmtx) + 6, label = sig1_high, hjust = "left") +
      # y-axis / sig2 annotations
      annotate("text", x = nrow(hypmtx) + 6, y = 0, label = sig2_low, hjust = "right", angle = 90) +
      annotate("text", x = nrow(hypmtx) + 6, y = nrow(hypmtx), label = sig2_high, hjust = "left", angle = 90) +
      coord_cartesian(clip = "off")

    if(waterfall){
      df_gg <- merged
      colnames(df_gg) <- c("id", "r1", "r2", "v1", "v2")
      df_gg$p1 <- ifelse(df_gg$v1 >= 0, T, F)
      df_gg$p2 <- ifelse(df_gg$v2 >= 0, T, F)
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
        filename = paste0(savename,"_heatmap.png"),
        plot = plt,
        width = wd, height = ht
      )
    }

    print(plt)
  }

  if(scatter){
    metricsct <- Rubrary::plot_scatter(
      df = merged,
      xval = paste0(metric1, ".1"), xlabel = sig1_name,
      yval = paste0(metric2, ".2"), ylabel = sig2_name,
      title = paste0("Metric Scatter Plot"),
      reverse = F
    ) +
      scale_x_continuous(breaks = scales::pretty_breaks(6)) +
      scale_y_continuous(breaks = scales::pretty_breaks(6))

    print(metricsct)

    ranksct <- Rubrary::plot_scatter(
      df = merged, rank = T, reverse = F,
      xval = paste0(metric1, ".1"), xlabel = sig1_name,
      yval = paste0(metric2, ".2"), ylabel = sig2_name,
      title = paste0("Rank Scatter Plot")
    ) +
      theme(line = element_blank(),
            axis.text = element_blank(),
            axis.title.x = element_text(size = rel(1.25), margin = margin(t = -(0.004 * nrow(merged)))),
            axis.title.y = element_text(size = rel(1.25), margin = margin(r = -(0.004 * nrow(merged))))) +
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

    print(ranksct)

    if(!is.null(savename)){
      ggsave(
        filename = paste0(savename,"_metricsct.png"),
        plot = metricsct,
        width = 4, height = 4
      )

      ggsave(
        filename = paste0(savename,"_ranksct.png"),
        plot = ranksct,
        width = wd, height = ht + 1
      )
    }
  }

  return(plt)
}
