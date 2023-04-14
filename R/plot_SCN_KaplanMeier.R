utils::globalVariables(c(
  "time", "estimate",
  "conf.low", "strata"
))

#' Plot SCN Kaplan Meier Survival Analysis
#'
#' @import ggplot2
#'
#' @param df dataframe; cols "OS" (overall survival censor), "OS.time" (overall survival time)
#' @param xlab string; x-axis label
#' @param title plot title
#' @param SCN SCN threshold to designate SCN-like status if not already in DF
#' @param savename string; filename to save as
#'
#' @return Kaplan-Meier overall survival analysis plot comparing SCN-like status
#' @export
#'
plot_SCN_KaplanMeier <- function(df, xlab = "Time", title = "Survival Analysis",
                                 SCN = NULL, savename = NULL) {
  Rubrary::use_pkg("survival", "ggsurvfit")

    if (!is.null(SCN)) {
      df$SCN.like <- "Non-SCN"
      df[df$Zscored_SCN_score > SCN, ]$SCN.like <- "SCN-like"
    } else {
      SCN <- 3
    }

  sf <- ggsurvfit::survfit2(survival::Surv(OS.time, OS) ~ SCN.like, data = df)
  pval <- ggsurvfit::survfit2_p(sf)
  km_plot <- sf %>%
    ggsurvfit::tidy_survfit() %>%
    ggplot(aes(
      x = time, y = estimate,
      ymin = conf.low, ymax = conf.low,
      linetype = strata
    )) +
    geom_step() +
    labs(
      x = xlab,
      y = "Overall Survival",
      title = title,
      subtitle = paste0("SCN-like: Z-scored SCN score > ", SCN)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    scale_x_continuous(breaks = seq(0, 8000, by = 2000), expand = c(0.02, 0)) +
    theme_classic() +
    scale_linetype_manual(
      values = c("solid", "longdash"),
      labels = c(
        paste0(
          "Non-SCN (",
          length(which(df$SCN.like == "Non-SCN")), ")"
        ),
        paste0(
          "SCN-like (",
          length(which(df$SCN.like == "SCN-like")), ")"
        )
      )
    ) +
    theme(
      legend.position = c(0.7, 0.8),
      legend.title = element_blank()
    ) +
    annotate("text", x = 4000, y = 0.7, label = paste0("Log-rank ", pval))

  if (!is.null(savename)) {
    ggsave(savename, km_plot, width = 8, height = 8)
  }

  message(survival::coxph(survival::Surv(OS.time, OS) ~ SCN.like, data = df))
  return(km_plot)
}
