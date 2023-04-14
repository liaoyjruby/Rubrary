
#' Prompt user to install package
#'
#' Helper function for "Rubrary::use_pkg"
#'
#' @param pkg string; package name (CRAN or Bioconductor)
#'
#' @return None
inst_pkg <- function(pkg){
  if (utils::menu(c("Yes", "No"), title = paste0("\nInstall package ", pkg, "?")) == "1") {
    BiocManager::install(pkg)
  } else { message(paste0("** ", pkg, " not installed; function may break!"))}
}

#' Check if package is installed, prompt if not
#'
#' Based on https://stackoverflow.com/a/44660688
#'
#' @param ... string; package(s) to check if installed
#'
#' @return None
#'
#' @examples use_pkg("ggplot2", "dplyr")
#' @export
use_pkg <- function(...){
  pkgs <- unlist(list(...))
  check <- invisible(unlist(lapply(pkgs, require,
                                   character.only = TRUE, quietly = TRUE)))
  need <- pkgs[check==FALSE]
  if(length(need) > 0){
    warning(
      paste0("Required packages not found: ", paste(need, collapse = ", ")),
      immediate. = TRUE
      )
    if (!require("BiocManager", quietly = TRUE)){
      inst_pkg("BiocManager")
    }
    invisible(lapply(need, inst_pkg))
  }
}

#' Reload package by installation path
#'
#' @param pkg string; name of package
#'
#' @return Reloads pkg
#' @export
reload_pkg <- function(pkg = "Rubrary"){
  devtools::reload(pkgload::inst(pkg))
}

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
  dfmerged <- dplyr::left_join(
    tibble::rownames_to_column(df1, var = "rn"),
    tibble::rownames_to_column(df2, var = "rn"), by = "rn")
  return(tibble::column_to_rownames(dfmerged, var = "rn"))
}

#' Full join by rownames
#'
#' Wrapper for dplyr::full_join to merge by rownames.
#'
#' @param df1 Numerical dataframe with rownames
#' @param df2 Numerical dataframe with rownames
#'
#' @return Full-joined dataframe by rownames with rownames
#' @seealso [dplyr::full_join()], [left_join_rownames()]
#' @export
full_join_rownames <- function(df1, df2){
  dfmerged <- dplyr::full_join(
    tibble::rownames_to_column(df1, var = "rn"),
    tibble::rownames_to_column(df2, var = "rn"), by = "rn")
  return(tibble::column_to_rownames(dfmerged, var = "rn"))
}

#' Combine data.frames by column, filling in missing rows.
#'
#' `cbinds` a list of dataframes filling missing rows with NA.
#' Based on https://stackoverflow.com/a/7962980
#'
#' @param ... input data frames to column bind together
#'
#' @return a single data frame
#' @export
#'
cbind.fill <- function(...) {
  trans <- lapply(list(...),t)
  trans_dfs <- lapply(trans, as.data.frame)
  return (data.frame(t(plyr::rbind.fill(trans_dfs))))
}


#' Like head but just the corner
#'
#' @param df dataframe
#' @param n integer; number of rows = columns to display
#'
#' @return Top left n x n submatrix of df
#' @export
#'
#' @examples
#' corner(diag(15))
corner <- function(df, n = 10){
  nr <- ifelse(n > nrow(df), nrow(df), n)
  nc <- ifelse(n > ncol(df), ncol(df), n)
  return(df[1:nr, 1:nc])
}
