utils::globalVariables(c(
  "..density.."
))

#' Plot distribution / normality check
#'
#' @param df dataframe
#' @param value string; colname of values
#' @param swtest logical; perform SW test for normality?
#' @param title string; plot title
#' @param xlab string; x-axis label
#' @param save logical; save as PNG?
#' @param savename string; file to save figure as
#'
#' @return Plot of distribution with median + mean indicated
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline labs theme_classic ggsave
#' @importFrom stats median shapiro.test
#' @export
#'
#' @examples
#' set.seed(13)
#' df_norm <- data.frame(Value = rnorm(50))
#' plot_distribution(df_norm, value = "Value", title = "Normal", swtest = TRUE)
#'
#' df_uni <- data.frame(Value = 1:100)
#' plot_distribution(df_uni, value = "Value", title = "Uniform", swtest = TRUE)
#'
plot_distribution <- function(df = NA, value, swtest = FALSE,
                              title = "Distribution", xlab = "Value",
                              save = FALSE, savename = "Distribution.png") {
  if (!is.character(value)) {
    df <- as.data.frame(as.numeric(value))
    colnames(df)[1] <- "Value"
    value <- "Value"
  }

  if (swtest) {
    # Shapiro-Wilk test of normality
    sw <- shapiro.test(df[, value])
    pval <- signif(sw$p.value, digits = 3)
    message(paste0(sw$method, ": p-value = ", pval))
    if (sw$p.value > 0.05) {
      msg_sw <- paste0("* p = ", pval, " > 0.05; normal")
    } else {
      msg_sw <- paste0("* p = ", pval, " <= 0.05; non-normal")
    }
    message(msg_sw)
  }

  plt <- ggplot(df, aes_string(x = value)) +
    geom_histogram(aes(y = ..density..), color = "black", fill = "white") +
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
