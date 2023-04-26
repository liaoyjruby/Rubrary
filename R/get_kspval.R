
#' Calculate Kolgomorov-Smirnov enrichment p-value
#'
#' If no group of interest specified, pulls first condition against others.
#'
#' @param df dataframe
#' @param value string; colname of values to rank by
#' @param group string; colname of categorizations
#' @param group2 string; group of interest
#'
#' @return Kolgomorov-Smirnov enrichment p-value
#' @export
#'
#' @examples
#' set.seed(13)
#' df = data.frame(
#'   group = c(rep("A", 50), rep("B",50)),
#'   values = c(rnorm(50, mean = 0), rnorm(50, mean = 2)))
#' get_kspval(df, value = "values", group = "group", group2 = "A")

get_kspval <- function(df, value, group, group2 = NULL){
  df <- df[order(df[, value]), ]
  df$rank <- 1:nrow(df)

  if(is.null(group2)) {
    group2 = unique(df[, group])[1]
  }

  ks_pval <- stats::ks.test(
    df[df[, group] == group2, "rank"],
    df[!(df[, group] == group2), "rank"]
  )$p.value

  return(ks_pval)
}
