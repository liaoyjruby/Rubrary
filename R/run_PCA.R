utils::globalVariables(c(
  "PC", "value", "variable"
))

#' Plot screeplot from prcomp
#'
#' @param obj_prcomp prcomp function output object
#' @param savename string; filepath to save PNG under
#' @param label logical; T for % values at points
#'
#' @return Screeplot
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_y_continuous xlab ylab labs theme_bw theme ggsave
#' @export
#'
plot_screeplot <- function(obj_prcomp, savename = NA, label = F){
  var_explained = obj_prcomp$sdev^2 / sum(obj_prcomp$sdev^2)
  df_scrplt <- data.frame(PC = colnames(obj_prcomp$rotation),
                          Var.Exp = var_explained[1:length(colnames(obj_prcomp$rotation))],
                          Cum.Var.Exp = cumsum(var_explained)[1:length(colnames(obj_prcomp$rotation))])

  df_scrplt$PC <- factor(df_scrplt$PC, levels = unique(df_scrplt$PC))
  df_scrplt <- reshape2::melt(df_scrplt, id.var = "PC")
  df_scrplt$value <- df_scrplt$value * 100
  df_scrplt$label <- paste0(round(df_scrplt$value, 1), "%")

  scrplt <- ggplot(data = df_scrplt, aes(x = PC, y = value, col = variable, group = variable, label = label)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    xlab("Principal Component") +
    ylab("Variance Explained (%)") +
    labs(title = paste0(basename(savename), " Screeplot")) +
    theme_bw() +
    theme(legend.title = element_blank()) +
    {if (label) ggrepel::geom_text_repel(size = 3, show.legend = F)}

  plot(scrplt)

  if(!is.na(savename)){
    ggsave(
      filename = paste0(savename,"_prcomp_screeplot.png"),
      plot = scrplt,
      width = 8, height = 6
    )
  }

  return(scrplt)
}

#' Run PCA (prcomp wrapper)
#'
#' Adapted from glab.library "PCA_from_file"
#'
#' @param df numeric dataframe; samples as columns, genes as rows
#' @param savename string; directory path and file name but no extension pls
#' @param summary logical; output summary info
#' @param center logical; indicate whether the variables should be shifted to be zero centered
#' @param scale logical; indicate whether the variables should be scaled to have unit variance before the analysis takes place
#' @param tol numerical; indicate the magnitude below which components should be omitted
#' @param screeplot logical; output + save screeplot?
#'
#' @return prcomp obj; output text files of PCA scores, loadings, sdev
#' @importFrom utils write.table
#'
#' @export
#'
run_PCA <- function(df, savename = "PCA", summary = FALSE,
                    center = TRUE, scale = TRUE, tol = 0.05, screeplot = TRUE){

  if (is.character(df)){
    dfpath <- df
    df <- read.delim(dfpath, row.names = 1)
  }

  # Filter out rows of all zero
  df <- df[rowSums((df[,-1]==0))<ncol(df[-1]),]

  # TRANSPOSE
  t.df=t(df)

  pca <- stats::prcomp(t.df, center=center, scale = scale, tol = tol)

  pca_scores <- pca$x
  pca_scores <- cbind("Score"=rownames(pca_scores),pca_scores)
  pca_loadings <- pca$rotation
  pca_loadings=cbind("Loading"=colnames(t.df),pca_loadings)
  pca_evalues=pca$sdev

  write.table(pca_scores,
              paste0(savename,"_prcomp_scores.txt"),
              sep='\t',row.names=FALSE,quote=FALSE)
  write.table(pca_loadings,
              paste0(savename,"_prcomp_loadings.txt"),
              sep='\t',row.names=FALSE,quote=FALSE)
  write.table(pca_evalues,
              paste0(savename,"_prcomp_sdev.txt"),
              sep='\t',row.names=FALSE,quote=FALSE)

  if (screeplot) {
    plot_screeplot(
      obj_prcomp = pca,
      savename = savename
    )
  }

  return(pca)
}
