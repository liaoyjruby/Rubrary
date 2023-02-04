utils::globalVariables(c(
  "..density.."
))

#' Shapiro-Wilk test of normality
#'
#' @param values numeric vector
#'
#' @return Logical SW normality test result; message
#' @export
#'
#' @examples
#' set.seed(13)
#' check_normal(rnorm(50))
#'
check_normal <- function(values){
  sw <- stats::shapiro.test(values)
  pval <- signif(sw$p.value, digits = 3)
  message(paste0(sw$method, ": p-value = ", pval))
  if (sw$p.value > 0.05) {
    msg_sw <- paste0("* p = ", pval, " > 0.05; normal")
    result <- TRUE
  } else {
    msg_sw <- paste0("* p = ", pval, " <= 0.05; non-normal")
    result <- FALSE
  }
  message(msg_sw)
  return(result)
}


#' Plot distribution
#'
#' @param df dataframe
#' @param value string; colname of values
#' @param check_normal logical; perform SW test for normality?
#' @param hist logical; include histogram?
#' @param title string; plot title
#' @param xlab string; x-axis label
#' @param save logical; save as PNG?
#' @param savename string; file to save figure as
#'
#' @return Plot of distribution with median + mean indicated
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline labs theme_classic ggsave
#' @importFrom stats median
#' @export
#'
#' @examples
#' set.seed(13)
#' df_norm <- data.frame(Value = rnorm(50))
#' plot_distribution(df_norm, value = "Value", title = "Normal", check_normal = TRUE)
#'
#' df_uni <- data.frame(Value = 1:100)
#' plot_distribution(df_uni, value = "Value", title = "Uniform", check_normal = TRUE)
#'
plot_distribution <- function(df = NA, value, check_normal = FALSE, hist = TRUE,
                              title = "Distribution", xlab = "Value",
                              save = FALSE, savename = "Distribution.png") {
  if (!is.character(value)) {
    df <- as.data.frame(as.numeric(value))
    colnames(df)[1] <- "Value"
    value <- "Value"
  }

  if (check_normal) {check_normal(df[,value])}

  plt <- ggplot(df, aes_string(x = value)) +
    {if(hist) geom_histogram(aes(y = ..density..), color = "black", fill = "white")} +
    geom_density(alpha = 0.2, fill = "red") +
    geom_vline(aes(xintercept = mean(df[, value])), color = "blue", linetype = "dashed") +
    geom_vline(aes(xintercept = median(df[, value])), color = "red", linetype = "dashed") +
    labs(
      title = title,
      subtitle = paste0("Median: ", round(median(df[, value]), 2), "; Mean: ", round(mean(df[, value]), 2)),
      x = xlab,
      y = "Density"
    ) +
    theme_classic()

  if (save) {
    ggsave(
      filename = savename,
      plot = plt,
      height = 5, width = 7
    )
  }

  plt
}
