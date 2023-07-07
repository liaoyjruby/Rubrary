
utils::globalVariables(c(
  "PC", "value", "variable",
  "Highlight", ".", "var"
))

#' Extract loadings from `prcomp` object
#'
#' Scales the PCA decomposition *eigenvectors* (`rotation` output of `prcomp`) by the square root of corresponding eigenvalues (`sdev` output of `prcomp`).
#'
#' @param obj_prcomp `prcomp` object
#'
#' @return Matrix of PCA loadings
#' @seealso `stats::prcomp`
#' @export
#'
#' @examples
#' data(iris)
#' PCA_iris <- Rubrary::run_PCA(t(iris[,c(1:4)]))
#' Rubrary::get_loadings(PCA_iris)
#'
get_loadings <- function(obj_prcomp){
  # `sdev` output from prcomp is sqrt(eigenvalues) of cor/cov mtx
  # diag(obj_prcomp$sdev)` is a square mtx w/ sqrt eigenvalues on diagonal
  # these are ACTUAL loadings, NOT eigenvectors
  rotmtx <- obj_prcomp$rotation
  loadings <- rotmtx %*% diag(obj_prcomp$sdev[1:ncol(rotmtx)]) %>%
    as.data.frame() %>%
    `colnames<-`(sub("V", "PC", colnames(.))) %>%
    as.matrix()
  return(loadings)
}

#' Run PCA (`prcomp` wrapper)
#'
#' Adapted from `glab.library::PCA_from_file`.
#'
#' In general, Z-score standardization (`center = T`; `scale = T`) before PCA is advised. For (transformed) gene expression data, genearlly, center but don't scale.
#'
#' `center = T`: PCA maximizes the sum-of-squared deviations *from the origin* in the first PC. Variance is only maximized if the data is pre-centered.
#'
#' `scale = T`: If one feature varies more than others, the feature will dominate resulting principal components. Scaling will also result in components in the same order of magnitude.
#'
#' Use either `tol` or `rank`, but not both.
#'
#' @import dplyr
#'
#' @param df (path to) numeric dataframe; samples as columns, genes/features as rows
#' @param summary logical; output summary info
#' @param center logical; indicate whether the variables should be shifted to be zero centered
#' @param scale logical; indicate whether the variables should be scaled to have unit variance
#' @param tol numeric; indicate the magnitude below which components should be omitted
#' @param rank integer; a number specifying the maximal rank, i.e., maximal number of principal components to be used
#' @param screeplot logical; output + save screeplot?
#' @param savename string; filepath (no ext.) to save PCA scores, loadings, sdev under
#'
#' @return prcomp obj
#' @seealso [stats::prcomp()], [Rubrary::plot_PCA()]
#' @export
#'
#' @examples
#' data(iris)
#' Rubrary::run_PCA(t(iris[,c(1:4)]))
#'
run_PCA <- function(df, savename = NULL, summary = FALSE,
                    center = TRUE, scale = FALSE, tol = 0.05, rank = NULL, screeplot = TRUE){

  if (is.character(df)){ df <- Rubrary::rread(df, row.names = 1) }
  t.df=t(df) # Transpose
  pca <- stats::prcomp(t.df, center=center, scale = scale, tol = tol, rank. = rank)
  if(summary){summary(pca)}

  pca_scores <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Scores")

  # pca_eigenvectors <- pca$rotation %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column(var = "Eigenvectors")

  pca_loadings <- Rubrary::get_loadings(pca) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Loadings")

  # `sdev` output from prcomp is sqrt(eigenvalues) of cor/cov mtx
  pca_sdev <- pca$sdev %>%
    as.data.frame()

  if(!is.null(savename)){
    Rubrary::rwrite(pca_scores,
                paste0(savename,"_prcomp_scores.txt"))
    # Rubrary::rwrite(pca_eigenvectors,
    #             paste0(savename,"_prcomp_eigenvectors.txt")
    Rubrary::rwrite(pca_loadings,
                paste0(savename,"_prcomp_loadings.txt"))
    Rubrary::rwrite(pca_sdev,
                paste0(savename,"_prcomp_sdev.txt"))
  }

  if (screeplot) {
    scrplt <- Rubrary::plot_screeplot(
      obj_prcomp = pca,
      savename = savename
    )
    print(scrplt)
  }

  return(pca)
}

#' Plot PCA with annotation
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param df_pca string or `prcomp` obj; (path to) PCA output
#' @param anno string or df; Annotation info for DF
#' @param PCx string; Component on x-axis
#' @param PCy string; Component on y-axis
#' @param type c("Score", "Loading")
#' @param label logical; T to label points
#' @param annoname string; Colname in `anno` matching point name
#' @param annolabel string; Colname in `anno` to label points by, defaults to `annoname`
#' @param annotype string; Colname in `anno` with info to color by
#' @param annotype2 string; Colname in `anno` with info to change shape by
#' @param ellipse logical; Draw `ggplot2::stat_ellipse` data ellipse w/ default params - this is NOT a confidence ellipse
#' @param ks_pval c("none", "caption", "grob"); Display ks-pvalue as caption or grob
#' @param title string; Plot title
#' @param subtitle string; Subtitle for plot
#' @param density logical; Show density plot along both axes; requires group annotations to be provided
#' @param highlight char vector; Specific points to shape differently & label
#' @param colors char vector; For discrete `annotype`, length should be number of unique `annotype`s. For continuous `annotype`, can either be length 2 where `colors[1]` is low and `colors[2]` is high or length 3 diverging colorscale where `colors[1]` = low, `colors[2]` = mid, `colors[3]` = high.
#' @param col_midpt numeric; For continuous `length(colors) == 3` `scale_color_gradientn` usage only
#' @param savename string; File path to save plot under
#' @param height numeric; Saved plot height
#' @param width numeric; Saved plot width
#'
#' @return PCA output plotted with annotation
#' @export
#' @examples
#' data(iris)
#' iris$Sample = rownames(iris)
#' PCA_iris <- Rubrary::run_PCA(t(iris[,c(1:4)]))
#' # Scores
#' Rubrary::plot_PCA(
#'   df_pca = PCA_iris,
#'   anno = iris[,c("Sample", "Species")],
#'   annoname = "Sample", annotype = "Species",
#'   title = "Iris PCA Scores by Species",
#'   subtitle = "Centered & scaled",
#'   ellipse = TRUE
#' )
#' # Loadings
#' Rubrary::plot_PCA(
#'   df_pca = PCA_iris,
#'   type = "Loadings",
#'   title = "Iris PCA Loadings",
#'   subtitle = "Centered & scaled",
#'   label = TRUE
#' )
#'
plot_PCA <- function(
    df_pca, anno = NULL, PCx = "PC1", PCy = "PC2", type = c("Scores", "Loadings"),
    annoname = "Sample", annolabel = annoname, label = FALSE,
    annotype = "Type", annotype2 = NULL, ellipse = FALSE,
    ks_pval = c("none", "caption", "grob"), highlight = NULL, colors = NULL, col_midpt = 0,
    title = NULL, subtitle = NULL, density = FALSE,
    savename = NULL, width = 8, height = 8) {
  type <- match.arg(type)
  ks_pval <- match.arg(ks_pval)
  annolabel <- ifelse(annolabel == annoname, type, annolabel)

  if (is.character(df_pca)) { # PCA results as path to txt
    dfpath <- df_pca
    df_pca <- Rubrary::rread(dfpath)
    sdevpath <- gsub("_[^_]+$", "_sdev.txt", dfpath)
    sdev <- Rubrary::rread(sdevpath) %>%
      mutate(var = .^2, #
             pve = round(var / sum(var) * 100, digits = 2)) %>%
      `rownames<-`(paste0("PC", rownames(.)))

    if (grepl("score", dfpath, ignore.case = T)) {
      type <- "Scores"
    } else if(grepl("loading", dfpath, ignore.case = T)){
      type <- "Loadings"
    }
  } else if (methods::is(df_pca, "prcomp")){ # prcomp object
    df_sc <- df_pca$x %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Scores")

    df_lo <- Rubrary::get_loadings(df_pca) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Loadings")

    if(type == "Loadings"){
      df <- df_lo
    } else { # Scores or biplot
      df <- df_sc
    }

    sdev <- df_pca$sdev %>%
      as.data.frame() %>%
      mutate(var = .^2,
             pve = round(var / sum(var) * 100, digits = 2)) %>%
      `rownames<-`(paste0("PC", rownames(.)))

  } else { # PCA results directly as dataframe
    df <- df_pca
    names(df)[1] <- type
  }

  if(exists("sdev")){
    PCxlab <- paste0(PCx, " (", sdev$pve[match(PCx, rownames(sdev))], "%)")
    PCylab <- paste0(PCy, " (", sdev$pve[match(PCy, rownames(sdev))], "%)")
  } else {
    PCxlab <- PCx
    PCylab <- PCy
  }

  if (is.character(anno)) {
    anno <- Rubrary::rread(anno)
  }

  if(!is.null(anno)){
    df <- df %>%
      left_join(., anno, by = stats::setNames(nm = type, annoname))
  }

  if (all(is.null(anno))) { # No coloring / annotation
    if(!is.null(colors)){
      col = colors[1]
    } else {
      col = "black"
    }

    plt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]])) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.25) +
      {if (type == "Scores" || type == "Biplot") geom_point(size = 2, color = col)} +
      {if (type == "Loadings") geom_segment(
        data = df, mapping = aes(x = 0, y = 0, xend = .data[[PCx]], yend = .data[[PCy]]),
        arrow = arrow(length = unit(0.025, "npc")), color = col)} +
      labs(title = title,
           x = PCxlab,
           y = PCylab) +
      {if (!is.null(subtitle)) labs(subtitle = subtitle)} +
      {if (label) ggrepel::geom_text_repel(label = df[, type], max.overlaps = Inf)} +
    theme_classic()
  } else {
    # KS pvalue calculation ----
    ks_cap <- waiver()
    if ((length(unique(df[, annotype])) == 2) && (ks_pval != "none")) {
      kpX <- Rubrary::get_kspval(df, PCx, annotype, unique(df[, annotype])[1])
      kpY <- Rubrary::get_kspval(df, PCy, annotype, unique(df[, annotype])[1])
      kpX_text <- paste0(
        sub("comp", "C", PCx), ": KS enrich. p = ",
        format.pval(kpX, digits = 3, eps = 0.001, nsmall = 3)
      )
      kpY_text <- paste0(
        sub("comp", "C", PCy), ": KS enrich. p = ",
        format.pval(kpY, digits = 3, eps = 0.001, nsmall = 3)
      )
      if (ks_pval == "grob"){
        grobX <- grid::grobTree(
          grid::textGrob(kpX_text,
            x = 0.95, y = 0.95, just = "right",
            gp = grid::gpar(col = "black", fontsize = 15, fontface = "italic")
          )
        )
        grobY <- grid::grobTree(
          grid::textGrob(kpY_text,
            x = 0.95, y = 0.9, just = "right",
            gp = grid::gpar(col = "black", fontsize = 15, fontface = "italic")
          )
        )
      }
      ks_cap <- ifelse(ks_pval == "caption", paste0(kpX_text, "; ", kpY_text), NULL)
    }

    # Manage alpha if many points
    alpha <- ifelse(nrow(df) > 10000, 0.7, 1)

    # Manage colors
    if(is.null(colors) && !is.numeric(df[, annotype])){
      Rubrary::use_pkg("scales")
      cols = scales::hue_pal()(length(unique(df[, annotype])))
    } else if (!is.null(colors) && colors[1] == "alpha"){
      Rubrary::use_pkg("Seurat") # Use Seurat implementation of pals instead
      cols = Seurat::DiscretePalette(n = length(unique(df[, annotype])), palette = "alphabet")
    } else {
      cols = colors
    }

    # Manage legend title + colors if numeric
    if (is.numeric(df[, annotype])) {
      guidetitle <- guide_colorbar(title = annotype)
      if(is.null(colors)){ cols <- c("blue", "red") } else { cols <- colors }
      if(length(cols) == 1){ cols <- c("gray80", cols) }
      if(length(cols) == 2){
        cols_numeric <- scale_color_gradient( # Numeric anno
          low = cols[1], high = cols[2], guide = "colourbar")
      } else if(length(cols) >= 3){
        cols_numeric <- scale_color_gradient2( # Numeric anno
          low = cols[1], mid = cols[2], high = cols[3],
          midpoint = col_midpt, guide = "colourbar")
      }
      density = F
      ellipse = F # NO ELLIPSE ALLOWED
    } else {
      guidetitle <- guide_legend(title = annotype)
    }

    # Manage highlights / add'l annos
    if (!all(is.null(highlight))) {
      label <- F
      df$Highlight <- ifelse(df$Score %in% highlight, "HL", "")
      ggplt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]], shape = "Highlight"))
    } else if (!is.null(annotype2)) {
      ggplt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]], shape = .data[[annotype2]]))
    } else {
      ggplt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]]))
    }

    plt <- ggplt +
      # X/Y axes
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.25) +
      # Loadings arrow plot
      {if (type == "Loadings") geom_segment(
        data = df, mapping = aes(x = 0, y = 0, xend = .data[[PCx]], yend = .data[[PCy]]),
        arrow = arrow(length = unit(0.025, "npc")))} +
      # Scores scatter plot
      {if (type == "Scores") geom_point(aes(color = .data[[annotype]]), size = 2, alpha = alpha)} +
      # Color
      {if (!is.numeric(df[, annotype])) scale_color_manual(values=cols)} + # Categorical anno
      {if (is.numeric(df[, annotype])) cols_numeric } +
      # Ellipse
      {if (ellipse) stat_ellipse(
        data = df, aes(color = .data[[annotype]]))} +
      labs(
        title = title,
        caption = ks_cap,
        x = PCxlab,
        y = PCylab) +
      {if (!is.null(subtitle)) labs(subtitle = subtitle)} +
      guides(color = guidetitle) + # Set title of legend
      theme_classic() +
      {if (label) ggrepel::geom_text_repel(
        label = df[, annolabel], max.overlaps = Inf, color = "black")} +
      {if (!all(is.null(highlight))) ggrepel::geom_text_repel(
        aes(label = ifelse(Highlight == "HL", df[, annolabel], "")),
        color = "black", max.overlaps = Inf, size = 3)}
  }

  if (!is.null(savename)) {
    ggsave(
      plot = plt,
      filename = savename,
      height = height,
      width = width
    )
  }

  if (density && !is.null(anno)) {
    Rubrary::use_pkg("ggExtra")

    if ((length(unique(df[, annotype])) == 2) && (ks_pval == "grob")) {
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
    plt <- mplt
  }

  return(plt)
}

#' Plot screeplot from prcomp
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param obj_prcomp prcomp function output object
#' @param savename string; filepath to save PNG under
#' @param label logical; T for % values at points
#' @param npcs integer; # of PCs to include in plot
#' @param cum_var_exp integer; c(1:100), cumulative variance explained threshold
#'
#' @return Screeplot as ggplot object
#' @export
#' @examples
#' data(iris)
#' PCA_iris <- Rubrary::run_PCA(t(iris[,c(1:4)]), screeplot = FALSE)
#' Rubrary::plot_screeplot(PCA_iris)
#'
plot_screeplot <- function(obj_prcomp, npcs = ncol(obj_prcomp$x), label = FALSE,
                           cum_var_exp = 80, savename = NULL){
  var_explained <- obj_prcomp$sdev^2 / sum(obj_prcomp$sdev^2)
  df_scrplt <- data.frame(
    PC = colnames(obj_prcomp$rotation),
    Var.Exp = var_explained[1:length(colnames(obj_prcomp$rotation))],
    Cum.Var.Exp = cumsum(var_explained)[1:length(colnames(obj_prcomp$rotation))]) %>%
    mutate(PC = as.numeric(sub("PC", "", PC)))

  message(paste0("** Cumulative var. exp. >= ", cum_var_exp,
                 "% at PC ",df_scrplt[df_scrplt$Cum.Var.Exp >= (cum_var_exp/100),]$PC[1],
                 " (", round(df_scrplt[df_scrplt$Cum.Var.Exp >= (cum_var_exp/100),]$Cum.Var.Exp[1] * 100, 1), "%)"))

  df_scrplt2 <- df_scrplt[1:npcs,] %>%
    tidyr::pivot_longer(!PC) %>%
    mutate(value = value * 100,
           label = paste0(round(value, 1), "%"))

  title_scr <- ifelse(is.null(savename), "PCA", basename(savename))

  scrplt <- ggplot(df_scrplt2, aes(x = PC, y = value, col = name, group = name, label = label)) +
    geom_hline(yintercept = cum_var_exp, color = "gray", linetype = "dashed") +
    geom_line() +
    geom_point() +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    xlab("Principal Component") +
    ylab("Variance Explained (%)") +
    labs(title = paste0(title_scr, " Screeplot")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          # axis.text.x = element_text(angle = 90)
    ) +
    {if (label) ggrepel::geom_text_repel(size = 3, show.legend = F)}

  if(!is.null(savename)){
    ggsave(
      filename = paste0(savename,"_prcomp_screeplot.png"),
      plot = scrplt,
      width = 8, height = 6
    )
  }

  return(scrplt)
}

#' Plot PCA matrix via `GGally::ggpairs`
#'
#' Only does "scores" PCA plots.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param df_pca string/`prcomp` obj; (path to) PCA output with sample names in first column
#' @param PCs num vector; list of numeric PCs to plot (ex. `c(1:3)` for first 3 PCs)
#' @param anno string/df; Annotation info for DF
#' @param annoname string; Colname in `anno` matching point name
#' @param annotype string; Colname in `anno` with info to color by
#' @param title string; Plot title
#' @param colors char vector; For discrete `annotype`, length should be number of unique `annotype`s.
#' @param savename string; File path to save plot under
#' @param height numeric; Saved plot height
#' @param width numeric; Saved plot width
#'
#' @return `ggmatrix` object with PCA scatter plot matrix
#' @export
#'
#' @examples
#' data(iris)
#' iris$Sample = rownames(iris)
#' plot_PCA_matrix(
#'   df_pca = Rubrary::run_PCA(t(iris[,c(1:4)]), screeplot = FALSE),
#'   PCs = c(1:3),
#'   anno = iris[,c("Sample", "Species")],
#'   annoname = "Sample",
#'   annotype = "Species",
#'   title = "Iris PCA"
#' )
plot_PCA_matrix <- function(
    df_pca, PCs = c(1:3), anno = NULL, annoname = "Sample", annotype = "Type",
    colors = NULL, title = NULL, savename = NULL, width = 8, height = 8) {

  # Load data ----
  if (is.character(df_pca)) { # PCA results as path to txt
    dfpath <- df_pca
    df_pca <- Rubrary::rread(dfpath)
    sdevpath <- gsub("_[^_]+$", "_sdev.txt", dfpath)
    sdev <- Rubrary::rread(sdevpath) %>%
      mutate(var = .^2, #
             pve = round(var / sum(var) * 100, digits = 2)) %>%
      `rownames<-`(paste0("PC", rownames(.)))
  } else if (methods::is(df_pca, "prcomp")){ # prcomp object
    df <- df_pca$x %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Scores")

    sdev <- df_pca$sdev %>%
      as.data.frame() %>%
      mutate(var = .^2,
             pve = round(var / sum(var) * 100, digits = 1)) %>%
      `rownames<-`(paste0("PC", rownames(.)))

  } else { # PCA results directly as dataframe, no sdev
    df <- df_pca
    names(df)[1] <- "Scores"
  }

  PCs <- paste0("PC", PCs)
  idx <- match(PCs, names(df))

  if(exists("sdev")){
    PClabs <- paste0(PCs, " (", sdev$pve[match(PCs, rownames(sdev))], "%)")
  } else { PClabs <- PCs }

  # Join annotation ----
  if (is.character(anno)) { anno <- Rubrary::rread(anno) }
  if(!is.null(anno)){
    df <- df %>% left_join(., anno, by = stats::setNames(nm = "Scores", annoname))
    axlab = "show"
    leg = c(length(idx), 1)
    legpos = "bottom"
  } else { #No anno
    df <- df %>% mutate(anno = "none")
    annotype = "anno"
    axlab = "internal"
    leg = NULL
    legpos = "none"
    if(is.null(colors)){ colors = "black" } else { colors = colors[1] }
  }

  # Plot mtx ----
  Rubrary::use_pkg("GGally", strict = TRUE)
  if(is.numeric(df[,annotype])){ # Continuous annotype

    guidetitle <- guide_colorbar(title = annotype)
    if(is.null(colors)){ cols <- c("blue", "red") } else { cols <- colors }
    if(length(cols) == 1){ cols <- c("gray80", cols) }

    num_sct <- function(data, mapping, ..., low = "blue", high = "red"){
      ggplot(data = data, mapping = mapping) +
        geom_point(...) +
        scale_color_gradient(low = low, high = high)
    }

    plt <- GGally::ggpairs(
      data = df, aes(color = .data[[annotype]]),
      columns = idx,
      title = title,
      upper = list(continuous = GGally::wrap(num_sct, low = cols[1], high = cols[2])),
      lower = list(continuous = GGally::wrap(num_sct, low = cols[1], high = cols[2])),
      # diag = list(continuous = GGally::wrap(GGally::ggally_diagAxis, color = "black")),
      axisLabels = "internal",
      columnLabels = PClabs,
      legend = leg,
      progress = FALSE
    ) +
      theme_bw() +
      theme(legend.position = legpos,
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  } else { # Discrete annotype
    guidetitle <- guide_legend(title = annotype)
    if(is.null(colors)){ colors <- scales::hue_pal()(length(unique(df[[annotype]])))}

    plt <- GGally::ggpairs(
      data = df, aes(color = .data[[annotype]], fill = .data[[annotype]]),
      columns = idx,
      title = title,
      columnLabels = PClabs,
      upper = list(continuous = "points"),
      diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5)),
      axisLabels = axlab,
      legend = leg,
      progress = FALSE
    ) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_classic() +
      theme(legend.position = legpos,
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
  }


  # Save ----
  if(!is.null(savename)){
    ggsave(
      plot = plt,
      filename = savename,
      width = width, height = height
    )
  }
  return(plt)
}

#' Plot PCA biplot
#'
#' For PCA SVD \eqn{X = USV^T}, uses standardized principal components for data points (\eqn{\bf{U}\sqrt{n - 1}}) and loadings (\eqn{\bf{VS}/\sqrt{n-1}}) and plots onto the same scale - a "proper" PCA biplot according to [this Stack Exchange thread](https://stats.stackexchange.com/questions/141085/positioning-the-arrows-on-a-pca-biplot/141531) citing the [Gabriel 1971 paper on PCA biplots](https://doi.org/10.2307/2334381).
#'
#' @param obj `prcomp` object
#' @param PCx string; Component on x-axis
#' @param PCy string; Component on y-axis
#' @param anno df; Annotation info for observations
#' @param annoname string; Colname in `anno` matching data points
#' @param annotype string; Colname in `anno` for desired coloring
#' @param label c("Both", "Loadings", "Scores", "None"); what points to label
#' @param colors char vector; Length should be number of unique `annotype`s
#' @param col_load string; Color for loading arrow segments
#' @param title string; Plot title
#' @param ellipse logical; Draw `ggplot2::stat_ellipse` data ellipse w/ default params - this is NOT a confidence ellipse
#' @param savename string; File path to save plot under
#' @param height numeric; Saved plot height
#' @param width numeric; Saved plot width
#'
#' @return Biplot as `ggplot` object
#' @export
#'
#' @examples
#' data(iris)
#' iris$Sample = rownames(iris)
#' PCA_iris <- Rubrary::run_PCA(t(iris[,c(1:4)]),
#'   center = TRUE, scale = FALSE, screeplot = FALSE)
#' Rubrary::plot_PCA_biplot(
#'   obj = PCA_iris,
#'   anno = iris[,c("Sample", "Species")],
#'   annoname = "Sample", annotype = "Species",
#'   label = "Loadings", ellipse = TRUE, title = "Iris PCA Biplot")
plot_PCA_biplot <- function(obj, PCx = "PC1", PCy = "PC2",
                            anno = NULL, annoname = "Sample", annotype = "Batch",
                            label = c("Both", "Loadings", "Scores", "None"),
                            colors = NULL, col_load = "firebrick", title = NULL,
                            ellipse = FALSE, savename = NULL, height = 8, width = 8){
  label <- match.arg(label)

  if(methods::is(obj, "prcomp")){

  }
  df_sc <- obj$x %>%
    as.data.frame() %>%
    mutate(across(everything(), ~./sd(.))) %>% # Standardized PC scores, div by std dev for unit variance
    tibble::rownames_to_column(var = "Scores")

  if(!is.null(anno)){
    df_sc <- df_sc %>%
      left_join(., anno, by = stats::setNames(nm = "Scores", annoname))
  }

  df_lo <- Rubrary::get_loadings(obj) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Loadings")

  sdev <- obj$sdev %>%
    as.data.frame() %>%
    mutate(var = .^2,
           pve = round(var / sum(var) * 100, digits = 2)) %>%
    `rownames<-`(paste0("PC", rownames(.)))

  # Variance explained
  PCxlab <- paste0(PCx, " (", sdev$pve[match(PCx, rownames(sdev))], "%)")
  PCylab <- paste0(PCy, " (", sdev$pve[match(PCy, rownames(sdev))], "%)")

  # Manage colors
  if(is.null(colors) && !is.numeric(df_sc[, annotype])){
    Rubrary::use_pkg("scales")
    cols = scales::hue_pal()(length(unique(df_sc[, annotype])))
  } else if (!is.null(colors) && colors == "alpha"){
    Rubrary::use_pkg("Seurat")
    cols = Seurat::DiscretePalette(n = length(unique(df_sc[, annotype])), palette = "alphabet2")
  } else {
    cols = colors
  }

  plt <- ggplot(df_sc, aes(x = .data[[PCx]], y = .data[[PCy]])) + # Scores df input
    # X/Y axes
    geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.25) +
    geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.25) +
    # Loadings
    geom_segment(
      data = df_lo, mapping = aes(x = 0, y = 0, xend = .data[[PCx]], yend = .data[[PCy]]),
      arrow = arrow(length = unit(0.025, "npc")), color = col_load) +
    # Scores
    {if (is.null(anno)) geom_point()} +
    {if (!is.null(anno)) geom_point(
      data = df_sc,
      mapping = aes(color = .data[[annotype]])
    )} +
    # Scores data ellipse
    {if (ellipse) stat_ellipse(
      data = df_sc, aes(color = .data[[annotype]]))} +
    # Labels
    labs(title = title,
         x = PCxlab,
         y = PCylab) +
    {if (label != "Loadings" && label != "None")
      ggrepel::geom_text_repel(label = df_sc[, "Scores"], max.overlaps = Inf)} + # Scores
    {if (label != "Scores" && label != "None")
      ggrepel::geom_label_repel( # Loadings
        data = df_lo,
        label = df_lo[, "Loadings"], max.overlaps = Inf, size = 2.5)} + # Loadings only
    theme_classic()

  if (!is.null(savename)) {
    ggsave(
      plot = plt,
      filename = savename,
      height = height,
      width = width
    )
  }

  return(plt)
}

#' Plot 3D PCA scores/loadings via `plotly`
#'
#' Image saving requires `kaleido` python package setup via `reticulate` R package. See `?plotly::save_image` for more details and setup instructions. Pass in `savename` with `html` file extension to avoid these errors.
#'
#' @import dplyr
#'
#' @param df_pca string/`prcomp` obj; (path to) PCA output with sample names in first column
#' @param PCs num vector; list of numeric PCs to plot (ex. `c(1:3)` for first 3 PCs)
#' @param type `c("Score", "Loading")`
#' @param anno string/df; Annotation info for `df_pca` with `annoname`, `annotype`, and `annolabel` columns
#' @param annoname string; Colname in `anno` matching point name
#' @param annotype string; Colname in `anno` with info to color by
#' @param annolabel string; Colname in `anno` to label points by, defaults to `annoname`
#' @param label logical; T to label points
#' @param colors char vector; For discrete `annotype`, length should be number of unique `annotype`s or length 2 for continuous `annotype`s where `colors[1]` represents low values and `colors[2]` represents high values.
#' @param title string; Plot title
#' @param savename string; File path to save plot under, `html` if not an image format
#' @param height numeric; Saved plot height if saving as image format
#' @param width numeric; Saved plot width
#' @param df_pca string or `prcomp` obj; (path to) PCA output
#' @param rotate logical; T to have the HTML widget automatically rotate when opened. Only applicable if `savename` is not `NULL`
#'
#' @return `plotly` object
#' @export
#'
#' @examples
#' data(iris)
#' iris$Sample = rownames(iris)
#' plot_PCA_3D(
#'   df_pca = Rubrary::run_PCA(t(iris[,c(1:4)]), screeplot = FALSE),
#'   PCs = c(1:3),
#'   type = "Scores",
#'   anno = iris[,c("Sample", "Species")],
#'   annoname = "Sample", annotype = "Species",
#'   title = "Iris PCA Scores 3D"
#' )
#'
plot_PCA_3D <- function(
    df_pca, PCs = c(1:3), type = c("Scores", "Loadings"),
    anno = NULL, annoname = "Sample", annotype = "Type",
    annolabel = annoname, label = FALSE, colors = NULL,
    title = NULL, savename = NULL, rotate = FALSE,
    width = 10, height = 10){
  type <- match.arg(type)
  annolabel <- ifelse(annolabel == annoname, type, annolabel)

  # Load data ----
  if(is.character(df_pca)) { # PCA results as path to txt
    dfpath <- df_pca
    df_pca <- Rubrary::rread(dfpath)
    sdevpath <- gsub("_[^_]+$", "_sdev.txt", dfpath)
    sdev <- Rubrary::rread(sdevpath) %>%
      mutate(var = .^2, #
             pve = round(var / sum(var) * 100, digits = 2)) %>%
      `rownames<-`(paste0("PC", rownames(.)))

    if (grepl("score", dfpath, ignore.case = T)) {
      type <- "Scores"
    } else if(grepl("loading", dfpath, ignore.case = T)){
      type <- "Loadings"
    }
  } else if (methods::is(df_pca, "prcomp")){ # prcomp object
    df_sc <- df_pca$x %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Scores")

    df_lo <- Rubrary::get_loadings(df_pca) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Loadings")

    if(type == "Loadings"){
      df <- df_lo
    } else { # Scores or biplot
      df <- df_sc
    }
    sdev <- df_pca$sdev %>%
      as.data.frame() %>%
      mutate(var = .^2,
             pve = round(var / sum(var) * 100, digits = 2)) %>%
      `rownames<-`(paste0("PC", rownames(.)))
  } else { # PCA results directly as dataframe
    df <- df_pca
    names(df)[1] <- type
  }

  PCs <- paste0("PC", PCs)

  if(exists("sdev")){
    PClabs <- paste0(PCs, " (", sdev$pve[match(PCs, rownames(sdev))], "%)")
  } else { PClabs <- PCs }

  # Join annotation ----
  if(is.character(anno)) { anno <- Rubrary::rread(anno) }
  if(!is.null(anno) && type == "Scores"){
    df <- df %>% left_join(., anno, by = stats::setNames(nm = "Scores", annoname))
    if(is.null(colors)){ colors <- scales::hue_pal()(length(unique(df[[annotype]])))}
    showleg = TRUE
  } else { #No anno
    df <- df %>% mutate(anno = "none")
    annotype = "anno"
    showleg = FALSE
    if(is.null(colors)){ colors = "black" }
  }

  # Plot ----
  Rubrary::use_pkg("plotly", strict = TRUE)
  if(is.numeric(df[,annotype])){ # Continuous legend
    fig <- plotly::plot_ly(
      type = "scatter3d", mode = "markers",
      x = ~df[,PCs[1]],
      y = ~df[,PCs[2]],
      z = ~df[,PCs[3]],
      marker = list(color = ~df[,annotype],
                    colorbar = list(title = annotype),
                    colors = colors),
      showlegend = showleg)
    } else {
      fig <- plotly::plot_ly(
        type = "scatter3d", mode = "markers",
        x = ~df[,PCs[1]],
        y = ~df[,PCs[2]],
        z = ~df[,PCs[3]],
        color = ~df[,annotype],
        colors = colors,
        showlegend = showleg)
  }


  fig <- fig %>%
    plotly::add_trace(
      marker = list(
        size = 7,
        line = list(
          color = "black",
          opacity = 0.5,
          width = 1)),
      showlegend = F
    )
  if(label){ fig <- fig %>%
    plotly::add_text(text = ~df[, annolabel], showlegend = F)}

  # Rotation? ----
    if(rotate){
      # @ismirsehregal on StackOverflow
      # https://stackoverflow.com/questions/71042818/plot-ly-rotate-animation-in-r
      Rubrary::use_pkg("htmlwidgets")
      fig <- fig %>%
        plotly::layout(
          title = title,
          legend = list(title=list(text=annotype)),
          scene = list(
            camera = list(
              eye = list(
                x = 1.25,
                y = 1.25,
                z = 1.25
              ),
              center = list(x = 0, y = 0, z = 0),
            xaxis = list(title = PClabs[1]),
            yaxis = list(title = PClabs[2]),
            zaxis = list(title = PClabs[3])
            ))) %>%
        htmlwidgets::onRender("
            function(el, x){
        var id = el.getAttribute('id');
        var gd = document.getElementById(id);
        Plotly.update(id).then(attach);
        function attach() {
          var cnt = 0;

          function run() {
            rotate('scene', Math.PI / 180); # speed
            requestAnimationFrame(run);
          }
          run();

          function rotate(id, angle) {
            var eye0 = gd.layout[id].camera.eye
            var rtz = xyz2rtz(eye0);
            rtz.t += angle;

            var eye1 = rtz2xyz(rtz);
            Plotly.relayout(gd, id + '.camera.eye', eye1)
          }

          function xyz2rtz(xyz) {
            return {
              r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
              t: Math.atan2(xyz.y, xyz.x),
              z: xyz.z
            };
          }

          function rtz2xyz(rtz) {
            return {
              x: rtz.r * Math.cos(rtz.t),
              y: rtz.r * Math.sin(rtz.t),
              z: rtz.z
            };
          }
        };
      }
          ")
    } else {
      fig <- fig %>%
        plotly::layout(
          scene = list(
            xaxis = list(title = PClabs[1]),
            yaxis = list(title = PClabs[2]),
            zaxis = list(title = PClabs[3])),
          title = title,
          legend = list(title = list(text = annotype)))
    }
  # Save ----
  if(!is.null(savename)){
    # Save as HTML ----
    htmlname <- paste0(tools::file_path_sans_ext(savename), ".html")
    Rubrary::use_pkg("htmlwidgets")
    htmlwidgets::saveWidget(
      plotly::partial_bundle(fig),
      file = htmlname,
      selfcontained = TRUE)
    utils::browseURL(htmlname)

    ## Save as image ----
    img_fmt <- c("png", "jpeg", "jpg", "webp", "svg", "pdf")
    if(tools::file_ext(savename) %in% img_fmt){
      message("** Image saving requires `kaleido` python package setup via `reticulate` R package.")
      message("** See `?plotly::save_image` for more details.")
      Rubrary::use_pkg("reticulate")
      if(requireNamespace("reticulate", quietly = TRUE)){
        plotly::save_image(
          p = fig, file = savename,
          scale = 4)
      }
    }
  }
  return(fig)
}

#' Apply varimax rotation to PCA scores
#'
#' [StackExchange reference](https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r)
#'
#' @param obj_prcomp `prcomp` object
#' @param ncomp integer; number of components to perform rotation with
#' @param normalize logical; T for Kaiser normalization: rows scaled to unit length before rotation, then scaled back afterwards
#' @param savename string; filepath (no ext.) to save results under
#'
#' @return `varimax` object with varimax rotation mtx `rotation`, varimax rotated `loadings`, standardized scores `std_scores`, and scores `scores`
#' @seealso [stats::varimax()]
#' @export
#'
#' @examples
#' data(iris)
#' iris$Sample = rownames(iris)
#' PCA_iris <- Rubrary::run_PCA(t(iris[,c(1:4)]))
#' PCA_iris_varimax <- Rubrary::rotate_varimax(PCA_iris)
#' head(PCA_iris_varimax$scores)
#' head(PCA_iris_varimax$loadings)
#' Rubrary::plot_PCA(
#'   df_pca = PCA_iris_varimax$scores,
#'   PCx = "V1", PCy = "V2",
#'   type = "Scores",
#'   anno = iris[,c("Sample", "Species")],
#'   annoname = "Sample", annotype = "Species",
#'   title = "Iris PCA Varimax"
#')
#'
rotate_varimax <- function(obj_prcomp, ncomp = 2, normalize = TRUE, savename = NULL){
  scores <- obj_prcomp$x[, 1:ncomp]
  loadings <- Rubrary::get_loadings(obj_prcomp)[, 1:ncomp]
  varimax_rotation <- stats::varimax(loadings, normalize = normalize)

  l <- varimax_rotation$loadings
  loadings_varimax <- data.frame(matrix(as.numeric(l), attributes(l)$dim, dimnames=attributes(l)$dimnames)) %>%
    rename_with(., ~ gsub("PC", "V", .x, fixed = TRUE)) %>%
    tibble::rownames_to_column(var = "Loadings")

  # Standardized scores
  std_scores_varimax <- scale(scores) %*% varimax_rotation$rotmat %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Scores")

  # "Unstandardize" scores by multiplying its values by sqrt(sum of squared loadings)
  ss_loadings <- colSums(loadings_varimax[, 2:(ncomp + 1)]^2)
  scores_varimax <- std_scores_varimax %>%
    tibble::column_to_rownames(var = "Scores") %>%
    as.matrix() %*% diag(sqrt(ss_loadings)) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Scores")

  if(!is.null(savename)){
    Rubrary::rwrite(
      loadings_varimax,
      file = paste0(savename, "_prcomp_loadings_varimax.txt"))
    Rubrary::rwrite(
      std_scores_varimax,
      file = paste0(savename, "_prcomp_std_scores_varimax.txt"))
    Rubrary::rwrite(
      scores_varimax,
      file = paste0(savename, "_prcomp_scores_varimax.txt"))
  }

  results_varimax <- list(
    rotation = varimax_rotation$rotmat,
    loadings = loadings_varimax,
    std_scores = std_scores_varimax,
    scores = scores_varimax
  )

  class(results_varimax) <- "varimax"

  return(results_varimax)
}

#' PCA prediction by projecting query/test samples onto PCA of reference/train samples
#'
#' Plot titles are automatically constructed based on parameters. Adding a `savedir` argument will result in intermediate plots and final projected scores to be output.
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom methods is
#'
#' @param train string/df; filepath to/df of reference/train samples data, genes/features as rownames, samples/observations as colnames
#' @param test string/df; filepath to/df of query/test samples data, genes/features as rownames, samples/observations as colnames
#' @param scale logical; T to scale variables to unit variance
#' @param varimax logical; T to perform varimax rotation on train + test
#' @param train_name string; descriptor for train samples
#' @param train_anno df; annotation info for train samples
#' @param train_annoname string; colname in `train_anno` matching point name
#' @param train_annotype string; colname in `train_anno` with info to color by
#' @param train_colors char vector; list of colors, length = # of unique `train_annotype`, set `NULL` for gray
#' @param test_name string; descriptor for test samples
#' @param test_anno df; annotation info for test samples
#' @param test_annoname string; colname in `test_anno` matching point name
#' @param test_annotype string; colname in `test_anno` with info to color by
#' @param test_colors char vector; list of colors, length = # of unique `test_annotype`
#' @param ellipse logical (vector); if length 2, `ellipse[1]` for train data ellipse, `ellipse[2]` for test data ellipse
#' @param label logical (vector); if length 2, `label[1]` for train data label, `label[2]` for test data label
#' @param flip NOT IMPLEMENTED - logical (vector); if length 2, `flip[1]` to flip x-axis values, `flip[2]` to flip y-axis values
#' @param savedir string; directory (+ prefix) to save output under; if directory, end string with "/"
#' @param height numeric; plot height
#' @param width numeric; plot width
#' @param fmt string; plot output format (ex. "png", "pdf")
#' @param test_only_plt logical; save scatter of projected test samples only
#' @param rank integer; maximal # of PCs to be used in `prcomp`; for compatibility w/ glab version
#'
#' @return ggplot object
#' @export
#' @seealso [Rubrary::run_PCA()], [Rubrary::plot_PCA()], [Rubrary::rotate_varimax()]
#' @examples
#' library(dplyr)
#' data(iris)
#' iris$Sample <- rownames(iris)
#' set.seed(13)
#' samp <- sample(nrow(iris), nrow(iris)*0.75) # Train data is 75% of iris
#' iris_train <- iris[samp,] %>%
#'   mutate(Batch = "Train",
#'          Species = paste0(Species, "_Train"))
#' iris_test <- iris[-samp,] %>%
#'   mutate(Batch = "Test",
#'          Species = paste0(Species, "_Test"))
#' Rubrary::predict_PCA(
#'   train = t(iris_train[,1:4]),
#'   test = t(iris_test[,1:4]),
#'   train_name = "Iris Train", train_anno = iris_train[,c("Sample", "Species")],
#'   train_annoname = "Sample", train_annotype = "Species",
#'   train_colors = c("pink", "red", "darkred"),
#'   test_name = "Iris Test", test_anno = iris_test[,c("Sample", "Species")],
#'   test_annoname = "Sample", test_annotype = "Species",
#'   test_colors = c("lightskyblue", "blue", "navy"),
#'   ellipse = c(FALSE, TRUE), label = FALSE
#' )
#'
predict_PCA <- function(
    train, test, scale = FALSE, varimax = FALSE,
    train_name = "Train", train_anno = NULL, train_annoname = NULL, train_annotype = NULL, train_colors = NULL,
    test_name = "Test", test_anno = NULL, test_annoname = NULL, test_annotype = NULL, test_colors = NULL,
    ellipse = FALSE, label = FALSE, flip = FALSE, savedir = NULL, height = 8, width = 8, fmt = "png",
    test_only_plt = FALSE, rank = 3){
  # Parameters ----
  PCx = "PC1"; PCy = "PC2" # Set as function arg?
  # Expand args to L2 if applicable
  if(length(ellipse) == 1){ ellipse = c(ellipse, ellipse)}
  ellipse_train <- ellipse[1]
  if(length(label) == 1){ label = c(label, label)}
  if(length(flip) == 1){ flip = c(flip, flip)}
  # Manage dataframes
  if(is(train, "data.frame") || is(train, "matrix")){ train_df <- train } else { train_df <- Rubrary::rread(train, row.names = 1)}
  if(is(test, "data.frame") || is(test, "matrix")){ test_df <- test } else { test_df <- Rubrary::rread(test, row.names = 1)}

  # Merge ----
  # Subset and order both train + test to common/intersect genes/features
  common_feats <- dplyr::intersect(rownames(train_df), rownames(test_df)) %>% sort()
  message(paste0("Common features: ", length(common_feats)))
  # Report unique features
  train_only_feats <- dplyr::setdiff(rownames(train_df), common_feats)
  test_only_feats <- dplyr::setdiff(rownames(test_df), common_feats)
  if(length(train_only_feats) != 0) {
    message(paste0("** ", train_name, " (train) unique features: ", length(train_only_feats),
                   " (", round((length(train_only_feats)/nrow(train_df)) * 100, 1),"%)"))
  }
  if(length(test_only_feats) != 0) {
    message(paste0("** ", test_name, " (test) unique features: ", length(test_only_feats),
                   " (", round((length(test_only_feats)/nrow(test_df)) * 100, 1),"%)"))
  }
  if(length(train_only_feats) != 0 || length(test_only_feats) != 0){
    message("All unique features excluded from train PCA & test projection!")
  }
  train_df <- train_df[common_feats,]
  test_df <- test_df[common_feats,]

  # Train PCA ----
  # Run train PCA
  train_pca <- Rubrary::run_PCA(
    df = train_df,
    center = TRUE, scale = scale,
    tol = NULL, rank = rank, screeplot = F)
  if(!varimax && flip[1]) { train_pca$x[,PCx] <- train_pca$x[,PCx] * -1}
  if(!varimax && flip[2]) { train_pca$x[,PCy] <- train_pca$x[,PCy] * -1}

  # Train annotation provided
  if(!is.null(train_anno)){
    if(is.character(train_anno)){ train_anno <- Rubrary::rread(train_anno) }
    # Filter anno to samples present
    train_anno <- train_anno %>%
      filter(!!as.symbol(train_annoname) %in% names(as.data.frame(train_df)))
    # No colors given or 1 color given
    if(is.null(train_colors) || length(train_colors) == 1){
      if(is.null(train_colors)){ train_colors <- "gray" }
      ellipse[1] <- FALSE # Use `ellipse_train` for train ellipse instead
      train_anno$ellipse <- train_anno[,train_annotype]
      train_anno[,train_annotype] <- train_name
    }
    train_pc_anno <- dplyr::left_join(
      x = tibble::rownames_to_column(as.data.frame(train_pca$x)[,c(PCx, PCy)], "Scores"),
      y = train_anno, by = stats::setNames(nm = "Scores", train_annoname))
  } else { # No train annotation
    label[1] <- FALSE
    ellipse[1] <- FALSE
    ellipse_train <- FALSE
    train_colors <- "gray"
    train_pc_anno <- data.frame(
      Samples = colnames(train_df),
      Type = train_name
    )
    train_annotype <- "Type"
  }

  train_plt <- Rubrary::plot_PCA(
    df_pca = train_pca,
    PCx = PCx, PCy = PCy,
    anno = train_anno,
    annoname = train_annoname,
    annotype = train_annotype,
    colors = train_colors,
    ellipse = ellipse[1],
    label = label[1],
    title = paste0("PCA - ", train_name)
  ) +
    {if(!is.null(train_anno) && !ellipse[1] && ellipse_train) stat_ellipse(
      data = train_pc_anno,
      mapping = aes(x = .data[[PCx]], y = .data[[PCy]], group = ellipse),
      color = train_colors)}

  if(!is.null(savedir)){
    ggsave(
      plot = train_plt,
      filename = paste0(savedir,"PCA_", gsub(" |\\/", "", train_name), ".", fmt),
      height = height, width = width
    )
  }
  # Rotate test data with train rotation mtx
  test_proj <- stats::predict(train_pca, newdata = t(test_df)) %>% as.data.frame()

  # Varimax ----
  if(varimax){
    ncomp <- 2 # Set as function arg?
    PCx = sub("PC", "V", PCx); PCy = sub("PC", "V", PCy)
    train_pca_vm <- Rubrary::rotate_varimax(train_pca, ncomp = ncomp, normalize = F)
    if(flip[1]) {train_pca_vm$scores[,PCx] <- train_pca_vm$scores[,PCx] * -1}
    if(flip[2]) {train_pca_vm$scores[,PCy] <- train_pca_vm$scores[,PCy] * -1}

    # Train annotation provided
    if(!is.null(train_anno)){
      train_pc_vm_anno <- dplyr::left_join(
        x = train_pca_vm$scores[,c("Scores", PCx, PCy)],
        y = train_anno, by = stats::setNames(nm = "Scores", train_annoname))
    } else {
      # No ellipse possible
      ellipse[1] <- FALSE
      ellipse_train <- FALSE
      train_colors <- "gray"
      train_pc_vm_anno <- data.frame(
        Samples = colnames(train_df),
        Type = train_name
      )
      train_annotype <- "Type"
    }

    train_plt_vm <- Rubrary::plot_PCA(
      df_pca = train_pca_vm$scores,
      PCx = PCx, PCy = PCy,
      anno = train_anno,
      annoname = train_annoname,
      annotype = train_annotype,
      colors = train_colors,
      ellipse = ellipse[1],
      label = label[1],
      title = paste0("PCA Varimax - ", train_name)
    ) +
      {if(!ellipse[1] && ellipse_train) stat_ellipse(
        data = train_pc_vm_anno,
        mapping = aes(x = .data[[PCx]], y = .data[[PCy]], group = ellipse),
        color = train_colors
      )}
    if(!is.null(savedir)){
      ggsave(
        plot = train_plt_vm,
        filename = paste0(savedir,"PCA_", gsub(" |\\/", "", train_name), "_varimax.", fmt),
        height = height, width = width
      )
    }
    train_plt <- train_plt_vm # Overwrite orig. plt w/ varimax version as base for projection
    # Test data varimax rotation
    test_proj_vm <- as.matrix(test_proj[,1:ncomp]) %*% train_pca_vm$rotation %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Scores")
    if(flip[1]) {test_proj_vm[,PCx] <- test_proj_vm[,PCx] * -1}
    if(flip[2]) {test_proj_vm[,PCy] <- test_proj_vm[,PCy] * -1}
    test_proj <- test_proj_vm
    train_pca <- train_pca_vm$scores
  } else {
    test_proj <- test_proj %>%
      tibble::rownames_to_column(var = "Scores")
  }

  # Test PCA ----
  if(!is.null(test_anno)){ # Test sample annotation present
    if(is.character(test_anno)){ test_anno <- Rubrary::rread(test_anno) } # Read if path
    # Filter to samples present
    test_anno <- test_anno %>%
      filter(!!as.symbol(test_annoname) %in% names(as.data.frame(test_df)))
    # Join
    test_proj_anno <- dplyr::left_join(
      test_proj[,c("Scores", PCx, PCy)], test_anno, by = stats::setNames(nm = "Scores", test_annoname))
    if(is.null(test_colors)){
      test_colors <- scales::hue_pal()(length(unique(test_proj_anno[,test_annotype])))
    }
  } else { # No test sample annotation
    # Set colors
    if(is.null(test_colors)){test_colors <- "red" } else { test_colors <- test_colors[1]}
    # Make fake merged test_anno
    test_annotype <- "Type"
    test_proj_anno <- test_proj %>%
      select(all_of(c("Scores", PCx, PCy))) %>%
      mutate(Type = test_name)
    label[2] <- FALSE
    ellipse[2] <- FALSE
  }
  colors <- stats::setNames(
    nm = c(as.character(sort(unique(train_pc_anno[,train_annotype]))),
           as.character(sort(unique(test_proj_anno[,test_annotype])))),
    c(train_colors, test_colors))

  vm_title <- ifelse(varimax, " Varimax", "")

  suppressMessages(
    test_plt <- train_plt +
      geom_point(data = test_proj_anno,
                 mapping = aes(x = .data[[PCx]], y = .data[[PCy]], color = .data[[test_annotype]])) +
      {if(ellipse[2]) stat_ellipse(data = test_proj_anno, aes(color = .data[[test_annotype]]))} +
      {if(label[2]) ggrepel::geom_text_repel(test_proj_anno, mapping = aes(x = .data[[PCx]], y = .data[[PCy]]),
                                             label = test_proj_anno$Scores, max.overlaps = Inf)} +
      scale_color_manual(values = colors, limits = names(colors), name = "Type") +
      guides(color = guide_legend(title = "Type")) + # Fix legend title
      {if(varimax) xlab(PCx)} + # Overwrite w/ "V1"
      {if(varimax) ylab(PCy)} + # Overwrite w/ "V2"
      labs(title = paste0("PCA", vm_title, " - ", test_name, " Projected on ", train_name))
  )

  # Save ----
  if(!is.null(savedir)){
    vm_name <- ifelse(varimax, "_varimax", "")
    # Save scores
    Rubrary::rwrite(
      x = test_proj,
      file = paste0(savedir, "PCA_proj_",
                    gsub(" |\\/", "", test_name), "_on_",
                    gsub(" |\\/", "", train_name), vm_name, ".txt"))
    # Save plot
    test_savename <- paste0(savedir, "PCA_proj_",
                            gsub(" |\\/", "", test_name), "_on_",
                            gsub(" |\\/", "", train_name), vm_name, ".", fmt)
    ggplot2::ggsave(
      filename = test_savename,
      plot = test_plt,
      height = height, width = width
    )
  }

  # Test only PCA ----
  if(test_only_plt){
    if(label[2]){
      label_test <- "Scores"
    } else {
      label_test <- NULL
    }
    if(!is.null(savedir)){
      test_only_savename <- paste0(
        savedir, "PCA_proj_", gsub(" |\\/", "", test_name), "_on_",
        gsub(" |\\/", "", train_name), vm_name, "_testonly.", fmt)
    } else {
      test_only_savename <- NULL
    }
    test_plt_only <- Rubrary::plot_scatter(
      df = test_proj_anno, xval = PCx, yval = PCy,
      label = label_test, group = test_annotype, colors = test_colors,
      cormethod = "none", guides = FALSE,
      title = paste0("PCA", vm_title, " - ", test_name, " Projected on ", train_name),
      subtitle = paste0(test_name, " samples only"),
      savename = test_only_savename,
      height = height, width = width
    )
  }
  return(test_plt)
}
