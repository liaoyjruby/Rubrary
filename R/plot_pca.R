utils::globalVariables(c(
  "Highlight"
))

#' Plot PCA with annotation
#'
#' @import ggplot2
#' @importFrom utils read.delim
#'
#' @param df string; Path to PC output
#' @param dfanno string or df; Annotation info for DF
#' @param PCx string; Component on x-axis
#' @param PCy string; Component on y-axis
#' @param PCtype c("Loading", "Score")
#' @param label logical; T to label points
#' @param annoname string; Colname in dfanno matching point name
#' @param annotype string; Colname in dfanno with info to color by
#' @param annotype2 string; Colname in dfanno with info to change shape by
#' @param title string; Plot title
#' @param subtitle string; Subtitle for plot
#' @param density logical; Show density plot along both axes
#' @param highlight char vector; Specific points to shape differently/labels
#' @param savename string; filepath to save plot under
#' @param width numeric; plot width
#' @param height numeric; plot width
#'
#' @return PCA output plotted with annotation
#'
#' @export
#'
plot_PCA <- function(df, dfanno = NULL, PCx = "PC1", PCy = "PC2", PCtype = "Score",
                     label = TRUE, annoname = "Sample", annotype = "Batch", annotype2 = NULL, highlight = NULL,
                     title = paste0("PCA ", PCtype, " Plot - ", annotype),
                     subtitle = NULL, density = FALSE, savename = NULL,
                     width = 8, height = 8) {
  if (is.character(df)) { # Read df from path
    dfpath <- df
    df <- read.delim(dfpath)
    sdevpath <- gsub("_[^_]+$", "_sdev.txt", dfpath)
    sdev <- read.delim(sdevpath)
    sdev$var <- unlist(sdev^2)
    sdev$pve <- unlist(round(sdev$var / sum(sdev$var) * 100, digits = 2))
    rownames(sdev) <- paste0("PC", seq(1, nrow(sdev)))

    if (grepl("score", dfpath, ignore.case = T)) {
      PCtype <- "Score"
    } else {
      PCtype <- "Loading"
    }
  }

  if(exists("sdev")){
    PCxlab <- paste0(PCx, " (", sdev$pve[match(PCx, rownames(sdev))], "%)")
    PCylab <- paste0(PCy, " (", sdev$pve[match(PCy, rownames(sdev))], "%)")
  } else {
    PCxlab <- PCx
    PCylab <- PCy
  }

  if (all(is.null(dfanno))) { # No coloring / annotation
    plt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]])) +
      geom_point(size = 2) +
      labs(
        title = title,
        x = PCxlab,
        y = PCylab
      ) +
      theme_classic() +
      {
        if (label) ggrepel::geom_text_repel(label = df[, PCtype])
      }
  } else {
    if (is.character(dfanno)) {
      dfanno <- read.delim(dfanno)
    }
    df[, 2:ncol(df)] <- sapply(df[, 2:ncol(df)], as.numeric)

    if (is.null(annotype2)) {
      df_merged <- dplyr::left_join(df, dfanno[, c(annoname, annotype)], by = structure(names = PCtype, .Data = annoname))
      df_merged <- df_merged[, c(PCtype, PCx, PCy, annotype)]
    } else {
      df_merged <- dplyr::left_join(df, dfanno[, c(annoname, annotype, annotype2)], by = structure(names = PCtype, .Data = annoname))
      df_merged <- df_merged[, c(PCtype, PCx, PCy, annotype, annotype2)]
    }
    # colnames(df_merged) <- c(PCx, PCy, "Type")

    # Two groups, calc KS p-value for both PCx and PCy
    if (length(unique(df_merged[, annotype])) == 2) {
      group1 <- unique(df_merged[, annotype])[1]
      # group2 <- unique(df_merged[,annotype])[2]

      # PCx ranking
      df_merged <- df_merged[order(df_merged[, PCx], decreasing = T), ]
      df_merged$rankX <- 1:nrow(df_merged)

      # PCy ranking
      df_merged <- df_merged[order(df_merged[, PCy], decreasing = T), ]
      df_merged$rankY <- 1:nrow(df_merged)

      kpX <- stats::ks.test( # for PCx
        x = df_merged[!df_merged[, annotype] == group1, "rankX"],
        y = df_merged[df_merged[, annotype] == group1, "rankX"]
      )$p.value

      kpY <- stats::ks.test( # for PCy
        x = df_merged[!df_merged[, annotype] == group1, "rankY"],
        y = df_merged[df_merged[, annotype] == group1, "rankY"]
      )$p.value

      grobX <- grid::grobTree(
        grid::textGrob(
          paste0(
            sub("comp", "C", PCx), ": KS enrich. p = ",
            format.pval(kpX, digits = 3, eps = 0.001, nsmall = 3)
          ),
          x = 0.95, y = 0.95, just = "right",
          gp = grid::gpar(col = "black", fontsize = 15, fontface = "italic")
        )
      )

      grobY <- grid::grobTree(
        grid::textGrob(
          paste0(
            sub("comp", "C", PCy), ": KS enrich. p = ",
            format.pval(kpY, digits = 3, eps = 0.001, nsmall = 3)
          ),
          x = 0.95, y = 0.9, just = "right",
          gp = grid::gpar(col = "black", fontsize = 15, fontface = "italic")
        )
      )

      df_merged <- df_merged[order(as.numeric(row.names(df_merged))), ]
    }

    if (nrow(df_merged) > 10000) {
      alpha <- 0.7
    } else {
      alpha <- 1
    }

    if (is.numeric(df_merged[, annotype])) {
      guidetitle <- guide_colorbar(title = annotype)
    } else {
      guidetitle <- guide_legend(title = annotype)
    }

    if (!all(is.null(highlight))) {
      label <- F
      df_merged$Highlight <- ifelse(df_merged$Score %in% highlight, "HL", "")
      ggplt <- ggplot(df_merged, aes(x = .data[[PCx]], y = .data[[PCy]], color = .data[[annotype]], shape = "Highlight"))
    } else if (!is.null(annotype2)) {
      ggplt <- ggplot(df_merged, aes(x = .data[[PCx]], y = .data[[PCy]], color = .data[[annotype]], shape = .data[[annotype2]]))
    } else {
      ggplt <- ggplot(df_merged, aes(x = .data[[PCx]], y = .data[[PCy]], color = .data[[annotype]]))
    }

    plt <- ggplt +
      geom_point(size = 2, alpha = alpha) +
      {
        if (is.numeric(df_merged[, annotype])) scale_color_gradient(low = "blue", high = "red", guide = "colourbar")
      } +
      labs(
        title = title,
        x = PCxlab,
        y = PCylab
      ) +
      {
        if (!is.null(subtitle)) labs(subtitle = subtitle)
      } +
      guides(color = guidetitle) + # Set title of legend
      theme_classic() +
      {
        if (label) ggrepel::geom_text_repel(label = df[, PCtype], max.overlaps = 30, color = "black")
      } + {
        if (!all(is.null(highlight))) ggrepel::geom_text_repel(aes(label = ifelse(Highlight == "HL", df[, PCtype], "")), color = "black", max.overlaps = 30, size = 3)
      }
  }

  if (!is.null(savename)) {
    ggsave(
      plot = plt,
      filename = savename,
      height = height,
      width = width
    )
  }

  if (density && !is.null(dfanno)) {
    if (length(unique(df_merged[, annotype])) == 2) {
      plt <- plt + annotation_custom(grobX) + annotation_custom(grobY)
    }
    mplt <- ggExtra::ggMarginal(plt,
      y = "Density", type = "density", margins = "both",
      size = 6, groupColour = T, groupFill = T
    )
    if (!is.null(savename)) {
      savename <- paste0(tools::file_path_sans_ext(savename), "_density.png")

      ggsave(
        plot = mplt,
        filename = savename,
        height = height,
        width = width
      )
    }
    # plt <- mplt
  }

  plt
}
