
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
  loadings <- obj_prcomp$rotation %*% diag(obj_prcomp$sdev) %>%
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
#' @importFrom utils read.delim
#'
#' @param df_pca string or prcomp obj; (path to) PCA output
#' @param anno string or df; Annotation info for DF
#' @param PCx string; Component on x-axis
#' @param PCy string; Component on y-axis
#' @param type c("Score", "Loading")
#' @param label logical; T to label points
#' @param annoname string; Colname in `anno` matching point name
#' @param annotype string; Colname in `anno` with info to color by
#' @param annotype2 string; Colname in `anno` with info to change shape by
#' @param title string; Plot title
#' @param subtitle string; Subtitle for plot
#' @param density logical; Show density plot along both axes
#' @param highlight char vector; Specific points to shape differently & label
#' @param colors char vector; length should be number of unique `annotype`s
#' @param savename string; filepath to save plot under
#' @param width numeric; plot width
#' @param height numeric; plot width
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
#'   title = "Iris PCA Scores by Species")
#' # Loadings
#' Rubrary::plot_PCA(df_pca = PCA_iris,
#'   type = "Loadings", title = "Iris PCA Loadings", label = TRUE)
#'
plot_PCA <- function(df_pca, anno = NULL, PCx = "PC1", PCy = "PC2", type = c("Scores", "Loadings"),
                     label = FALSE, annoname = "Sample", annotype = "Batch", annotype2 = NULL,
                     highlight = NULL, colors = NULL, title = NULL, subtitle = NULL, density = FALSE,
                     savename = NULL, width = 8, height = 8) {
  type <- match.arg(type)

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
    } else {
      type <- "Loadings"
    }
  } else if (methods::is(df_pca, "prcomp")){ # prcomp object
    if(type == "Scores"){
      df <- df_pca$x %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Scores")
    } else {
      df <- Rubrary::get_loadings(df_pca) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "Loadings")
    }

    sdev <- df_pca$sdev %>%
      as.data.frame() %>%
      mutate(var = .^2,
             pve = round(var / sum(var) * 100, digits = 2)) %>%
      `rownames<-`(paste0("PC", rownames(.)))

  } else {
    df <- df_pca
  }

  if(exists("sdev")){
    PCxlab <- paste0(PCx, " (", sdev$pve[match(PCx, rownames(sdev))], "%)")
    PCylab <- paste0(PCy, " (", sdev$pve[match(PCy, rownames(sdev))], "%)")
  } else {
    PCxlab <- PCx
    PCylab <- PCy
  }

  if (all(is.null(anno))) { # No coloring / annotation
    plt <- ggplot(df, aes(x = .data[[PCx]], y = .data[[PCy]])) +
      geom_point(size = 2) +
      labs(title = title,
           x = PCxlab,
           y = PCylab) +
      theme_classic() +
      {if (label) ggrepel::geom_text_repel(label = df[, type])}
  } else {
    if (is.character(anno)) {
      anno <- read.delim(anno)
    }
    df[, 2:ncol(df)] <- sapply(df[, 2:ncol(df)], as.numeric)

    if (is.null(annotype2)) {
      df_merged <- dplyr::left_join(df, anno[, c(annoname, annotype)],
                                    by = structure(names = type, .Data = annoname))
      df_merged <- df_merged[, c(type, PCx, PCy, annotype)]
    } else {
      df_merged <- dplyr::left_join(df, anno[, c(annoname, annotype, annotype2)],
                                    by = structure(names = type, .Data = annoname))
      df_merged <- df_merged[, c(type, PCx, PCy, annotype, annotype2)]
    }

    # Two groups, calc KS p-value for both PCx and PCy
    if (length(unique(df_merged[, annotype])) == 2) {
      kpX <- Rubrary::get_kspval(df_merged, PCx, annotype,unique(df_merged[, annotype])[1])
      kpY <- Rubrary::get_kspval(df_merged, PCy, annotype, unique(df_merged[, annotype])[1])

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

    # Manage alpha if many points
    alpha <- ifelse(nrow(df_merged) > 10000, 0.7, 1)

    # Manage colors
    if(is.null(colors)){
      Rubrary::use_pkg("scales")
      cols = scales::hue_pal()(length(unique(df_merged[, annotype])))
    } else if (colors == "alpha"){
      Rubrary::use_pkg("pals")
      cols = unname(pals::alphabet2(n = length(unique(df_merged[, annotype]))))
    } else {
      cols = colors
    }

    # Manage legend title + colors if numeric
    if (is.numeric(df_merged[, annotype])) {
      guidetitle <- guide_colorbar(title = annotype)
      if(is.null(colors)){cols <- c("blue", "red")}
    } else {
      guidetitle <- guide_legend(title = annotype)
    }

    # Manage highlights / add'l annos
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
      {if (!is.numeric(df_merged[, annotype])) scale_fill_manual(values=cols)} +
      {if (is.numeric(df_merged[, annotype])) scale_color_gradient(
        low = cols[1], high = , guide = "colourbar")} +
      labs(
        title = title,
        x = PCxlab,
        y = PCylab) +
      {if (!is.null(subtitle)) labs(subtitle = subtitle)} +
      guides(color = guidetitle) + # Set title of legend
      theme_classic() +
      {if (label) ggrepel::geom_text_repel(label = df[, type], max.overlaps = 30, color = "black")} +
      {if (!all(is.null(highlight))) ggrepel::geom_text_repel(
        aes(label = ifelse(Highlight == "HL", df[, type], "")), color = "black", max.overlaps = 30, size = 3)}
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

#' Apply varimax rotation to PCA scores
#'
#' [StackExchange reference](https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r)
#'
#' @param obj_prcomp `prcomp` object
#' @param ncomp integer; number of components to perform rotation with
#' @param normalize logical; T for Kaiser normalization: rows scaled to unit length before rotation, then scaled back afterwards
#' @param savename string; filepath (no ext.) to save results under
#'
#' @return List with varimax rotated loadings matrix and varimax rotated scores matrices
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
  scores_varimax <- scale(scores) %*% varimax_rotation$rotmat # Scores must be standardized / scaled
  loadings_varimax <- as.matrix(unlist(varimax_rotation$loadings))
  rownames(loadings_varimax) <- rownames(loadings)

  if(!is.null(savename)){
    write.table(
      tibble::rownames_to_column(scores, var = "Scores"),
      file = paste0(savename, "_prcomp_scores_varimax.txt"),
      quote = F, sep = "\t", row.names = F
    )
    write.table(
      tibble::rownames_to_column(varimax_rotation$loadings, var = "Loadings"),
      file = paste0(savename, "_prcomp_loadings_varimax.txt"),
      quote = F, sep = "\t", row.names = F
    )
  }

  results_rotated = list(
    loadings = loadings_varimax,
    scores = scores_varimax
  )

  return(results_rotated)
}
