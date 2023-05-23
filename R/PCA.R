
utils::globalVariables(c(
  "PC", "value", "variable",
  "Highlight", ".", "var",
  "kp_grob"
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
#' In general, Z-score standardization (`center = T`; `scale = T`) before PCA is advised.
#'
#' `center = T`: PCA maximizes the sum-of-squared deviations *from the origin* in the first PC. Variance is only maximized if the data is pre-centered.
#'
#' `scale = T`: If one feature varies more than others, the feature will dominate resulting principal components. Scaling will also result in components in the same order of magnitude.
#'
#' @import dplyr
#' @importFrom utils write.table
#'
#' @param df (path to) numeric dataframe; samples as columns, genes/features as rows
#' @param summary logical; output summary info
#' @param center logical; indicate whether the variables should be shifted to be zero centered
#' @param scale logical; indicate whether the variables should be scaled to have unit variance
#' @param tol numerical; indicate the magnitude below which components should be omitted
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
                    center = TRUE, scale = TRUE, tol = 0.05, screeplot = TRUE){

  if (is.character(df)){
    dfpath <- df
    df <- read.delim(dfpath, row.names = 1)
  }

  t.df=t(df) # Transpose

  pca <- stats::prcomp(t.df, center=center, scale = scale, tol = tol)
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
  pca_sdev <- pca$sdev

  if(!is.null(savename)){
    write.table(pca_scores,
                paste0(savename,"_prcomp_scores.txt"),
                sep='\t',row.names=FALSE,quote=FALSE)
    # write.table(pca_eigenvectors,
    #             paste0(savename,"_prcomp_eigenvectors.txt"),
    #             sep='\t',row.names=FALSE,quote=FALSE)
    write.table(pca_loadings,
                paste0(savename,"_prcomp_loadings.txt"),
                sep='\t',row.names=FALSE,quote=FALSE)
    write.table(pca_sdev,
                paste0(savename,"_prcomp_sdev.txt"),
                sep='\t',row.names=FALSE,quote=FALSE)
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
#' @importFrom utils read.delim
#'
#' @param df_pca string or `prcomp` obj; (path to) PCA output
#' @param anno string or df; Annotation info for DF
#' @param PCx string; Component on x-axis
#' @param PCy string; Component on y-axis
#' @param type c("Score", "Loading")
#' @param label logical; T to label points
#' @param annoname string; Colname in `anno` matching point name
#' @param annotype string; Colname in `anno` with info to color by
#' @param annotype2 string; Colname in `anno` with info to change shape by
#' @param ellipse logical; Draw `ggplot2::stat_ellipse` data ellipse w/ default params - this is NOT a confidence ellipse
#' @param ks_grob logical; Display ks-pvalue as grob instead of caption
#' @param title string; Plot title
#' @param subtitle string; Subtitle for plot
#' @param density logical; Show density plot along both axes
#' @param highlight char vector; Specific points to shape differently & label
#' @param colors char vector; Length should be number of unique `annotype`s
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
#' Rubrary::plot_PCA(df_pca = PCA_iris,
#'   anno = iris[,c("Sample", "Species")],
#'   annoname = "Sample", annotype = "Species",
#'   title = "Iris PCA Scores by Species",
#'   ellipse = TRUE)
#' # Loadings
#' Rubrary::plot_PCA(df_pca = PCA_iris,
#'   type = "Loadings", title = "Iris PCA Loadings", label = TRUE)
#'
plot_PCA <- function(df_pca, anno = NULL, PCx = "PC1", PCy = "PC2", type = c("Scores", "Loadings"),
                     label = FALSE, annoname = "Sample", annotype = "Batch", annotype2 = NULL,
                     ellipse = FALSE, ks_grob = FALSE, highlight = NULL, colors = NULL,
                     title = NULL, subtitle = NULL, density = FALSE,
                     savename = NULL, width = 8, height = 8) {
  type <- match.arg(type)
  lab_type <- type

  if (is.character(df_pca)) { # PCA results as path to txt
    dfpath <- df_pca
    df_pca <- read.delim(dfpath)
    sdevpath <- gsub("_[^_]+$", "_sdev.txt", dfpath)
    sdev <- read.delim(sdevpath) %>%
      as.data.frame() %>%
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
    anno <- read.delim(anno)
  }

  if(!is.null(anno)){
    df <- df %>%
      left_join(., anno, by = stats::setNames(nm = type, annoname))
  }

  if (all(is.null(anno))) { # No coloring / annotation
    plt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]])) +
      geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.25) +
      {if (type == "Scores" || type == "Biplot") geom_point(size = 2)} +
      {if (type == "Loadings") geom_segment(
        data = df, mapping = aes(x = 0, y = 0, xend = .data[[PCx]], yend = .data[[PCy]]),
        arrow = arrow(length = unit(0.025, "npc")))} +
      labs(title = title,
           x = PCxlab,
           y = PCylab) +
      {if (label) ggrepel::geom_text_repel(label = df[, lab_type])} +
    theme_classic()
  } else {
    # Two groups, calc KS p-value for both PCx and PCy
    ks_cap <- waiver()
    if (length(unique(df[, annotype])) == 2) {
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
      if (ks_grob){
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
      ks_cap <- ifelse(kp_grob, waiver(), paste0(kpX_text, "; ", kpY_text))
    }

    # Manage alpha if many points
    alpha <- ifelse(nrow(df) > 10000, 0.7, 1)

    # Manage colors
    if(is.null(colors) && !is.numeric(df[, annotype])){
      Rubrary::use_pkg("scales")
      cols = scales::hue_pal()(length(unique(df[, annotype])))
    } else if (!is.null(colors) && colors == "alpha"){
      Rubrary::use_pkg("pals")
      cols = unname(pals::alphabet2(n = length(unique(df[, annotype]))))
    } else {
      cols = colors
    }

    # Manage legend title + colors if numeric
    if (is.numeric(df[, annotype])) {
      guidetitle <- guide_colorbar(title = annotype)
      if(is.null(colors)){cols <- c("blue", "red")}
      density = F
      ellipse = F # NO ELLIPSE ALLOWED
    } else {
      guidetitle <- guide_legend(title = annotype)
    }

    # Manage highlights / add'l annos
    if (!all(is.null(highlight))) {
      label <- F
      df$Highlight <- ifelse(df$Score %in% highlight, "HL", "")
      ggplt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]], color = .data[[annotype]], shape = "Highlight"))
    } else if (!is.null(annotype2)) {
      ggplt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]], color = .data[[annotype]], shape = .data[[annotype2]]))
    } else {
      ggplt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]], color = .data[[annotype]]))
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
      {if (type == "Scores") geom_point(size = 2, alpha = alpha)} +
      # Color
      {if (!is.numeric(df[, annotype])) scale_color_manual(values=cols)} + # Categorical anno
      {if (is.numeric(df[, annotype])) scale_color_gradient( # Numeric anno
        low = cols[1], high = cols[2], guide = "colourbar")} +
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
      {if (label) ggrepel::geom_text_repel(label = df[, type], max.overlaps = 30, color = "black")} +
      {if (!all(is.null(highlight))) ggrepel::geom_text_repel(
        aes(label = ifelse(Highlight == "HL", df[, type], "")),
        color = "black", max.overlaps = 30, size = 3)}
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
    if (length(unique(df[, annotype])) == 2 && ks_grob) {
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
  var_explained = obj_prcomp$sdev^2 / sum(obj_prcomp$sdev^2)
  df_scrplt <- data.frame(PC = colnames(obj_prcomp$rotation),
                          Var.Exp = var_explained[1:length(colnames(obj_prcomp$rotation))],
                          Cum.Var.Exp = cumsum(var_explained)[1:length(colnames(obj_prcomp$rotation))])

  df_scrplt$PC <- sub("PC", "", df_scrplt$PC)
  message(paste0("** Cumulative var. exp. >= ", cum_var_exp,
                 "% at PC ",df_scrplt[df_scrplt$Cum.Var.Exp >= (cum_var_exp/100),]$PC[1],
                 " (", round(df_scrplt[df_scrplt$Cum.Var.Exp >= (cum_var_exp/100),]$Cum.Var.Exp[1] * 100, 1), "%)"))

  df_scrplt <- df_scrplt[1:npcs,]

  df_scrplt$PC <- factor(df_scrplt$PC, levels = unique(df_scrplt$PC))
  df_scrplt <- reshape2::melt(df_scrplt, id.var = "PC")
  df_scrplt$value <- df_scrplt$value * 100
  df_scrplt$label <- paste0(round(df_scrplt$value, 1), "%")

  title_scr <- ifelse(is.null(savename), "PCA", basename(savename))

  scrplt <- ggplot(data = df_scrplt, aes(x = PC, y = value, col = variable, group = variable, label = label)) +
    geom_hline(yintercept = cum_var_exp, color = "red") +
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

#' Plot "proper" PCA biplot
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
    Rubrary::use_pkg("pals")
    cols = unname(pals::alphabet2(n = length(unique(df_sc[, annotype]))))
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

#' Apply varimax rotation to PCA scores
#'
#' [StackExchange reference](https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r)
#'
#' @param obj_prcomp `prcomp` object
#' @param ncomp integer; number of components to perform rotation with
#' @param normalize logical; T for Kaiser normalization: rows scaled to unit length before rotation, then scaled back afterwards
#' @param savename string; filepath (no ext.) to save results under
#'
#' @return `varimax` object with varimax rotated `loadings`, standardized scores `std_scores`, and scores `scores`
#' @seealso [stats::varimax()]
#' @export
#'
#' @examples
#' data(iris)
#' PCA_iris <- Rubrary::run_PCA(t(iris[,c(1:4)]))
#' PCA_iris_varimax <- Rubrary::rotate_varimax(PCA_iris)
#' head(PCA_iris_varimax$scores)
#' head(PCA_iris_varimax$loadings)
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
    write.table(
      tibble::rownames_to_column(loadings_varimax$loadings, var = "Loadings"),
      file = paste0(savename, "_prcomp_loadings_varimax.txt"),
      quote = F, sep = "\t", row.names = F
    )
    write.table(
      std_scores_varimax,
      file = paste0(savename, "_prcomp_std_scores_varimax.txt"),
      quote = F, sep = "\t", row.names = F
    )
    write.table(
      scores_varimax,
      file = paste0(savename, "_prcomp_scores_varimax.txt"),
      quote = F, sep = "\t", row.names = F
    )
  }

  results_varimax <- list(
    loadings = loadings_varimax,
    std_scores = std_scores_varimax,
    scores = scores_varimax
  )

  class(results_varimax) <- "varimax"

  return(results_varimax)
}
