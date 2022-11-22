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
#'
#' @return output text files of PCA scores, loadings, sdev
#' @importFrom utils write.table
#'
#' @export
#'
run_PCA <- function(df, savename = "PCA", summary = FALSE,
                    center = TRUE, scale = TRUE, tol = 0.05){

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

  if (summary){
    print(summary(pca))
    stats::screeplot(pca)
    # pca
  }
}
