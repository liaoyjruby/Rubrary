
#' Wrapper for `data.table::fwrite` with tab separated values (TSV) default parameters
#'
#' @param x dataframe; `fwrite` input
#' @param file string; filepath to save output under, with extension
#' @param sep string; separator between columns
#' @param quote string; T to wrap fields in quotes
#' @param row.names logical; T to write rownames to file
#'
#' @return Dataframe written out to `savename` file
#' @export
rwrite <- function(x, file, sep = "\t", quote = FALSE, row.names = FALSE){
  data.table::fwrite(
    x = x,
    file = file,
    sep = sep,
    quote = quote,
    row.names = row.names
  )
}

#' Wrapper for `data.table::fread` implementing `row.names` functionality
#'
#' Default file reading function in Rubrary functions because data is often in large tables and `read.delim` isn't optimized for high dimensionality. If file extension is `xlsx`, will attempt to use `openxlsx::read.xlsx` instead. Setting `row.names != 0` will result in output being a `data.frame`.
#'
#' @import dplyr
#'
#' @param input string; `fread` input
#' @param row.names integer; column to use as row names, 0 for none
#' @param make.names logical; `TRUE` to ensure column names are unique and valid (replicating `utils::read.table`)
#' @param to_df logical; T to convert output to data.frame
#'
#' @return Table read in from input
#' @export
#'
#' @examples
#' glab_Beltran_2016 <- "https://raw.githubusercontent.com/graeberlab-ucla/glab.library/master/vignettes/PCA_tutorial/Beltran_2016_rsem_genes_upper_norm_counts_coding_log2.txt"
#' df <- Rubrary::rread(glab_Beltran_2016, row.names = 1, to_df = FALSE)
#' Rubrary::corner(df, 5)
rread <- function(input, row.names = 0, make.names = TRUE, to_df = TRUE){
  if(tools::file_ext(input) == "xlsx"){
    Rubrary::use_pkg("openxlsx")
    df <- openxlsx::read.xlsx(input)
  } else {
    df <- data.table::fread(input)
  }
  if(make.names){ colnames(df) <- make.names(colnames(df), unique = TRUE, allow_ = TRUE)}
  if(row.names != 0){ df <- tibble::column_to_rownames(df, var = names(df)[row.names])}
  if(to_df){df <- as.data.frame(df)}
  return(df)
}

#' Prompt user to install package
#'
#' Helper function for "Rubrary::use_pkg"
#' @keywords internal
#' @param pkg string; package name (CRAN or Bioconductor)
#' @param strict logical; strict requirement on package install?
#'
#' @return None
inst_pkg <- function(pkg, strict = FALSE){
  if (utils::menu(c("Yes", "No"), title = paste0("\nInstall package ", pkg, "?")) == "1") {
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      BiocManager::install(pkg)
    } else {
      utils::install.packages(pkg)
    }
  } else {
    if(strict){
      stop(paste0("** ", pkg, " not installed; stopping"))
    } else {
      message(paste0("** ", pkg, " not installed; function may not work as intended"))
    }
  }
}

#' Check if package is installed, prompt if not
#'
#' Based on https://stackoverflow.com/a/44660688
#'
#' @param ... string; package(s) to check if installed
#' @param strict logical; strict requirements on package install?
#'
#' @return None
#'
#' @examples use_pkg("ggplot2", "dplyr")
#' @export
use_pkg <- function(..., strict = FALSE){
  pkgs <- unlist(list(...))
  check <- invisible(unlist(lapply(pkgs, requireNamespace, quietly = TRUE)))
  need <- pkgs[check==FALSE]
  if(length(need) > 0){
    warning(
      paste0("Required packages not found: ", paste(need, collapse = ", ")),
      immediate. = TRUE
      )
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      inst_pkg("BiocManager")
    }
    invisible(lapply(need, inst_pkg, strict))
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

#' Split text into multiple lines
#'
#' @param text string
#' @param chars integer; rough # of characters per line
#' @param lines integer; desired # of lines
#'
#' @return `text` split into lines if longer than `split_nchar`
#' @export
#'
#' @examples
#' lipsum <- "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat."
#' writeLines(split_line(lipsum, chars = 50))
#' writeLines(split_line(lipsum, lines = 3))
#' writeLines(split_line("Short line of text!"))
#' writeLines(split_line("Short line of text!", lines = 2))
split_line <- function(text, chars = 40, lines = NULL){
  if(nchar(text) > chars || !is.null(lines)){
    words <- unlist(strsplit(text, " "))
    n_lines <- ifelse(
      is.null(lines), ceiling(nchar(text)/chars), lines)
    idx <- floor(seq(from = 1, to = length(words), length.out = n_lines + 1))
    idx <- sort(c(idx, idx[-1] + 1))
    lines <- c()
    for(i in 1:(n_lines*2)){
      if(i %% 2 == 1){ # Odd indices
        # message(paste0("Indices: ", i, " -> ", i+1))
        line <- paste(words[idx[i]:idx[i+1]], collapse = " ")
        # message(paste0("** Line: ", line))
        lines <- c(lines,line)
      }
    }
    text <- paste(lines, collapse = "\n")
  }
  return(text)
}
