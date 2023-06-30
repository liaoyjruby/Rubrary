
#' Calculate (signed) Kolgomorov-Smirnov enrichment p-value in dataframe
#'
#' If no group of interest specified, pulls first condition against others.
#'
#' Signed ES score adapted from [franapoli/signed-ks-test](https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R) which is based on [GSEA enrichment score calculation](https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/#enrichment-score-es).
#' Positive `ES` indicates that by rank, the `goi` are "left" of / smaller than the rest of the values, while negative `ES` indicates the `goi` is "right" of / greater than the rest of values. Use `viz = TRUE` to visualize the values and resulting KS p-value.
#'
#' @import dplyr
#'
#' @param df df; has columns `value` and `group`, `NA` values are filtered out
#' @param value string; colname of values to rank by
#' @param group string; colname of categorizations
#' @param goi group of interest value within `df[,group]`
#' @param rep0 numeric; value to replace `0` results with
#' @param signed logical; `TRUE` to return GSEA-style enrichment score as second element of list: `ks_pval[1] = pval`, `ks_pval[1] = ES`
#' @param viz logical; `TRUE` to output density plot of values
#'
#' @return Kolgomorov-Smirnov enrichment p-value & enrichment score
#' @export
#'
#' @examples
#' set.seed(13)
#' df = data.frame(
#'   group = c(rep("A", 50), rep("B",50)),
#'   values = c(rnorm(50, mean = 0), rnorm(50, mean = 2)))
#' Rubrary::get_kspval(df, value = "values", group = "group", goi = "A", signed = TRUE)

get_kspval <- function(df, value, group, goi = NULL, rep0 = .Machine$double.xmin, signed = FALSE, viz = FALSE){
  df <- df %>%
    filter(!is.na(!!sym(value))) %>%
    arrange(!!sym(value)) %>%
    mutate(rank = 1:nrow(.))

  if(is.null(goi)) {goi = unique(df[, group])[1]}

  pval <- stats::ks.test(
    df[df[, group] == goi, "rank"],
    df[!(df[, group] == goi), "rank"]
  )$p.value
  # Lowest possible value in base R - technically should be `.Machine$double.xmin`
  if(pval == 0) { pval <- rep0 }
  if(signed){ # GSEA-style enrichment value
    ks_pval <- list(pval = pval)
    # https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R
    tx <- df[df[, group] == goi, "rank"] # Ranks of `goi`
    n_tx <- length(tx) # Number of `goi`
    # bg <- which(!GSEA_pw_rnk$term)
    bg <- df[df[, group] != goi, "rank"] # Ranks of non-`goi`
    n_bg <- length(bg) # Number of non-`goi`
    zx <- cumsum(ifelse(order(c(tx, bg)) <= n_tx, 1/n_tx, -1/n_bg)) # Increment per pw
    ES <- zx[which.max(abs(zx))] # Enrichment score is max
    ks_pval$ES <- ES
  } else {
    ks_pval <- pval
    ES <- NULL
  }

  if(viz){
    ES_str <- ifelse(!is.null(ES), paste0("; ES = ", round(ES, 2)), NULL)
    plt <- Rubrary::plot_density(df = df, value = value, group = group, pval = F) +
      ggplot2::labs(
        title = paste0("GOI: ", goi),
        subtitle = paste0("KS enrich. p-value = ", signif(pval, digits = 4), ES_str)
      )
    print(plt)
  }
  return(ks_pval)
}

#' Shapiro-Wilk test of normality
#'
#' Uses `stats::shapiro.test` to perform test of normality and outputs interpretation into console.
#'
#' The null hypothesis of the S-W test is that the sample comes from a normally distributed population. If the p-value is less than the chosen alpha level, the null hypothesis is rejected, indicating that the data has non-normal distribution.
#'
#' @param values numeric vector
#' @param alpha numeric; p-value threshold for significance
#'
#' @return Logical SW normality test result; message
#' @seealso [Rubrary::plot_distribution()]
#' @export
#'
#' @examples
#' set.seed(13)
#' check_normal(rnorm(100))
#'
check_normal <- function(values, alpha = 0.05){
  if(length(values) > 5000){
    message("** Sampling 5000 values w/o replacement...")
    values <- sample(values, 5000)
  }
  sw <- stats::shapiro.test(values)
  pval <- signif(sw$p.value, digits = 3)
  message(paste0(sw$method, ": p-value = ", pval))
  result <- sw$p.value > alpha
  if (result) {
    msg_sw <- paste0("** p = ", pval, " > ", alpha, "; normal")
  } else {
    msg_sw <- paste0("** p = ", pval, " <= ", alpha, "; non-normal")
  }
  message(msg_sw)
  return(result)
}
