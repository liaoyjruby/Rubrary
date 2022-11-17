#' Left join by rownames
#'
#' Wrapper for dplyr::left_join to merge by rownames. Both dataframes should have matching rownames
#'
#' @param df1 Numerical dataframe with rownames
#' @param df2 Numerical dataframe with rownames
#'
#' @return Left-joined dataframe by rownames with rownames
#' @export
left_join_rownames <- function(df1, df2){
  dfmerged <- dplyr::left_join(tibble::rownames_to_column(df1, var = "rn"), tibble::rownames_to_column(df2, var = "rn"), by = "rn")
  rownames(dfmerged) <- dfmerged$rn
  dfmerged$rn <- NULL
  return(dfmerged)
}
