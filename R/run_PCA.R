#' Run PCA (prcomp wrapper)
#'
#' Adapted from glab.library "PCA_from_file"
#'
#' @param df (path to) numeric dataframe; samples as columns, genes as rows
#' @param savename string; filepath (no ext.) to save PCA scores, loadings, sdev under
#' @param summary logical; output summary info
#' @param center logical; indicate whether the variables should be shifted to be zero centered
#' @param scale logical; indicate whether the variables should be scaled to have unit variance before the analysis takes place
#' @param tol numerical; indicate the magnitude below which components should be omitted
#' @param screeplot logical; output + save screeplot?
#'
#' @return prcomp obj
#' @importFrom utils write.table
#'
#' @export
#'
run_PCA <- function(df, savename = NULL, summary = FALSE,
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

  if(!is.null(savename)){
    write.table(pca_scores,
                paste0(savename,"_prcomp_scores.txt"),
                sep='\t',row.names=FALSE,quote=FALSE)
    write.table(pca_loadings,
                paste0(savename,"_prcomp_loadings.txt"),
                sep='\t',row.names=FALSE,quote=FALSE)
    write.table(pca_evalues,
                paste0(savename,"_prcomp_sdev.txt"),
                sep='\t',row.names=FALSE,quote=FALSE)
  }

  if (screeplot) {
    Rubrary::plot_screeplot(
      obj_prcomp = pca,
      savename = savename
    )
  }

  return(pca)
}
