utils::globalVariables(c("density", "x"))

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


#' Plot distribution of numeric vector
#'
#' Ignores NA values for mean & median calculation.
#'
#' @import ggplot2
#'
#' @param values numeric vector; values to check distribution for
#' @param check_normal logical; perform SW test for normality?
#' @param hist logical; include histogram?
#' @param title string; plot title
#' @param xlab string; x-axis label
#' @param savename string; filepath to save figure under
#' @param height numeric; plot height
#' @param width numeric; plot width
#'
#' @return Plot of distribution with median + mean indicated
#' @seealso [Rubrary::check_normal()]
#' @export
#'
#' @examples
#' set.seed(13)
#' vals_normal <- rnorm(100)
#' plot_distribution(values = vals_normal, title = "Normal", hist = TRUE, check_normal = TRUE)
#'
#' vals_sequential <- c(1:100)
#' plot_distribution(values = vals_sequential, title = "Sequential", check_normal = TRUE)
#'
plot_distribution <- function(values, check_normal = FALSE, hist = FALSE,
                              title = "Distribution", xlab = "Value",
                              savename = NULL, height = 5, width = 7) {

  if (check_normal && (length(unique(values)) != 1)) {Rubrary::check_normal(values)}

  mea <- mean(as.numeric(values), na.rm = T)
  med <- stats::median(as.numeric(values), na.rm = T)

  plt <- ggplot(data.frame(x = values), aes(x = x)) +
    {if(hist) geom_histogram(aes(y = after_stat(density)), color = "black", fill = "white")} +
    geom_density(alpha = 0.2, fill = "red") +
    geom_vline(aes(xintercept = mea), color = "blue", linetype = "dashed") +
    geom_vline(aes(xintercept = med), color = "red", linetype = "dashed") +
    labs(
      title = title,
      subtitle = paste0("Median: ", round(med, 2), "; Mean: ", round(mea, 2)),
      x = xlab,
      y = "Density"
    ) +
    theme_classic()

  if (!is.null(savename)) {
    ggsave(
      filename = savename,
      plot = plt,
      height = height, width = width
    )
  }

  return(plt)
}
