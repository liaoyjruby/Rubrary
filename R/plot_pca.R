utils::globalVariables(c(
  "Highlight"
))

#' Plot PCA with annotation
#'
#' @param df string; Path to PC output
#' @param dfanno string or df; Annotation info for DF
#' @param PCx string; Component on x-axis
#' @param PCy string; Component on y-axis
#' @param PCtype c("Loading", "Score")
#' @param label logical; T to label points
#' @param annoname string; Colname in dfanno matching point name
#' @param annotype string; Colname in dfanno with annotation
#' @param title string; Plot title
#' @param marginal logical; Show density plot along both axes
#' @param save logical; T to save to png
#' @param savename string; Output plot name
#' @param highlight char vector; Specific points to shape differently/labels
#' @param subtitle string; Subtitle for plot
#'
#' @return PCA output plotted with annotation
#'
#' @importFrom ggplot2 ggplot aes aes_string geom_point labs theme_bw ggsave annotation_custom guide_legend guides guide_colorbar scale_color_gradient
#' @importFrom utils read.delim
#' @export
#'
plot_PCA <- function(df, dfanno = NA, PCx = "PC1", PCy = "PC2", PCtype = "Score",
                     label = T, annoname = "Sample", annotype = "Batch", highlight = NA,
                     title = paste0("PCA ",PCtype ," Plot - ", annotype),
                     subtitle = NA, marginal = F, save = F, savename = "PCA_plot.png"){

  if (is.character(df)){
    dfpath <- df
    df <- read.delim(dfpath)
    sdevpath <- gsub("_[^_]+$", "_sdev.txt", dfpath)
    sdev <- read.delim(sdevpath)
    sdev$var = unlist(sdev^2)
    sdev$pve = unlist(round(sdev$var/sum(sdev$var) * 100, digits = 2))
    rownames(sdev) = paste0("PC",seq(1,nrow(sdev)))

    if(grepl("score", dfpath, ignore.case = T)){
      PCtype = "Score"
    } else {
      PCtype = "Loading"
    }
  }

  if (all(is.na(dfanno))){ # No coloring / annotation
    plt <- ggplot(df, aes_string(x = PCx, y = PCy)) +
      geom_point(size = 2) +
      labs(title = title,
           x = paste0(PCx," (", sdev$pve[match(PCx, rownames(sdev))], "%)"),
           y = paste0(PCy," (", sdev$pve[match(PCy, rownames(sdev))], "%)")) +
      theme_classic() +
      {if(label) ggrepel::geom_text_repel(label = df[,PCtype])}
  } else {
    if (is.character(dfanno)){
      dfanno <- read.delim(dfanno)
    }
    df[,2:ncol(df)] <- sapply(df[,2:ncol(df)], as.numeric)

    df_merged <- dplyr::left_join(df, dfanno[,c(annoname, annotype)], by = structure(names = PCtype, .Data = annoname))
    df_merged <- df_merged[,c(PCtype, PCx, PCy, annotype)]
    # colnames(df_merged) <- c(PCx, PCy, "Type")

    # Two groups, calc KS p-value for both PCx and PCy
    if(length(unique(df_merged[,annotype])) == 2){
      group1 <- unique(df_merged[,annotype])[1]
      # group2 <- unique(df_merged[,annotype])[2]

      # PCx ranking
      df_merged <- df_merged[order(df_merged[,PCx], decreasing = T),]
      df_merged$rankX <- 1:nrow(df_merged)

      # PCy ranking
      df_merged <- df_merged[order(df_merged[,PCy], decreasing = T),]
      df_merged$rankY <- 1:nrow(df_merged)

      kpX <- stats::ks.test( # for PCx
        x = df_merged[!df_merged[,annotype] == group1, "rankX"],
        y = df_merged[df_merged[,annotype] == group1, "rankX"]
      )$p.value

      kpY <- stats::ks.test( # for PCy
        x = df_merged[!df_merged[,annotype] == group1, "rankY"],
        y = df_merged[df_merged[,annotype] == group1, "rankY"]
      )$p.value

      grobX <- grid::grobTree(
        grid::textGrob(paste0(sub("comp", "C", PCx), ": KS enrich. p = ",
                              format.pval(kpX, digits = 3, eps = 0.001, nsmall = 3)),
                       x = 0.95, y = 0.95, just = "right",
                       gp=grid::gpar(col="black", fontsize=15, fontface="italic")))

      grobY <- grid::grobTree(
        grid::textGrob(paste0(sub("comp", "C", PCy), ": KS enrich. p = ",
                              format.pval(kpY, digits = 3, eps = 0.001, nsmall = 3)),
                       x = 0.95, y = 0.9, just = "right",
                       gp=grid::gpar(col="black", fontsize=15, fontface="italic")))

      df_merged <- df_merged[order(as.numeric(row.names(df_merged))), ]
    }

    if (nrow(df_merged) > 10000){
      alpha = 0.7
    } else {
      alpha = 1
    }

    if(is.numeric(df_merged[,annotype])){
      guidetitle = guide_colorbar(title=annotype)
    } else {
      guidetitle = guide_legend(title=annotype)
    }

    if(!all(is.na(highlight))){
      label = F
      df_merged$Highlight <- ifelse(df_merged$Score %in% highlight, "HL", "")
      ggplt <- ggplot(df_merged, aes_string(x = PCx, y = PCy, color = annotype, shape = "Highlight"))

    } else {
      ggplt <- ggplot(df_merged, aes_string(x = PCx, y = PCy, color = annotype))
    }

    plt <- ggplt +
      geom_point(size = 2, alpha = alpha) +
      {if(is.numeric(df_merged[,annotype])) scale_color_gradient(low = "blue", high = "red", guide = "colourbar")} +
      labs(title = title,
           x = paste0(PCx," (", sdev$pve[match(PCx, rownames(sdev))], "%)"),
           y = paste0(PCy," (", sdev$pve[match(PCy, rownames(sdev))], "%)")) +
      {if(!is.na(subtitle)) labs(subtitle = subtitle)} +
      guides(color=guidetitle,
             shape = "none") + # Set title of legend
      theme_classic() +
      {if(label) ggrepel::geom_text_repel(label = df[,PCtype], max.overlaps = 30, color = "black")} +
      {if(!all(is.na(highlight))) ggrepel::geom_text_repel(aes(label=ifelse(Highlight == "HL",df[,PCtype],"")), color = "black", max.overlaps = 30, size = 3)}

  }

  if (save) {
    ggsave(
      plot = plt,
      filename = savename,
      height = 8,
      width = 8
    )
  }

  if (marginal && !is.na(dfanno)){
    if (length(unique(df_merged[,annotype])) == 2) {
      plt <- plt + annotation_custom(grobX) + annotation_custom(grobY)
    }
    mplt <- ggExtra::ggMarginal(plt,
                                y = "Density", type = "density", margins = "both",
                                size = 6, groupColour = T, groupFill = T)
    savename <- paste0(tools::file_path_sans_ext(savename), "_marginal.png")
    if (save){
      ggsave(
        plot = mplt,
        filename = savename,
        height = 8,
        width = 8
      )
    }
    plt <- mplt
  }

  plt
}
