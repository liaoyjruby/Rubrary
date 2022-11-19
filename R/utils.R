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

#' Full join by rownames
#'
#' Wrapper for dplyr::full_join to merge by rownames.
#'
#' @param df1 Numerical dataframe with rownames
#' @param df2 Numerical dataframe with rownames
#'
#' @return Full-joined dataframe by rownames with rownames
#' @export
full_join_rownames <- function(df1, df2){
  dfmerged <- dplyr::full_join(tibble::rownames_to_column(df1, var = "rn"), tibble::rownames_to_column(df2, var = "rn"), by = "rn")
  rownames(dfmerged) <- dfmerged$rn
  dfmerged$rn <- NULL
  return(dfmerged)
}

#' Combine data.frames by column, filling in missing rows.
#'
#' `cbinds` a list of dataframes filling missing rows with NA.
#'
#' @param ... input data frames to column bind together
#'
#' @return a single data frame
#' @export
#'
cbind.fill <- function(...) { # https://stackoverflow.com/a/7962980
  transpoted <- lapply(list(...),t)
  transpoted_dataframe <- lapply(transpoted, as.data.frame)
  return (data.frame(t(plyr::rbind.fill(transpoted_dataframe))))
}
