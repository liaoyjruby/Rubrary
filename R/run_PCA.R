#' Run PCA (prcomp wrapper)
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
#' @param savename string; filepath (no ext.) to save PCA scores, loadings, sdev under
#' @param summary logical; output summary info
#' @param center logical; indicate whether the variables should be shifted to be zero centered
#' @param scale logical; indicate whether the variables should be scaled to have unit variance
#' @param tol numerical; indicate the magnitude below which components should be omitted
#' @param screeplot logical; output + save screeplot?
#'
#' @return prcomp obj
#' @seealso [stats::prcomp()]
#'
#' @export
#'
run_PCA <- function(df, savename = NULL, summary = FALSE,
                    center = TRUE, scale = TRUE, tol = 0.05, screeplot = TRUE){

  if (is.character(df)){
    dfpath <- df
    df <- read.delim(dfpath, row.names = 1)
  }

  df <- df[rowSums((df[,-1]==0))<ncol(df[-1]),] # Filter out rows of all zero
  t.df=t(df) # Transpose

  pca <- stats::prcomp(t.df, center=center, scale = scale, tol = tol)

  pca_scores <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Score")

  pca_loadings <- pca$rotation %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Loading")

  pca_evalues <- pca$sdev

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
    scrplt <- Rubrary::plot_screeplot(
      obj_prcomp = pca,
      savename = savename
    )
    print(scrplt)
  }

  return(pca)
}
