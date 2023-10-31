utils::globalVariables(c("density", "x"))

#' Plot distribution of numeric vector
#'
#' Ignores NA values for mean & median calculation.
#'
#' @import ggplot2
#'
#' @param values numeric vector; values to check distribution for
#' @param check_normal logical; perform SW test for normality?
#' @param ctr_measure logical; include central measures (mean + median)?
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
plot_distribution <- function(values, check_normal = FALSE,
                              ctr_measure = TRUE, hist = FALSE,
                              title = "Distribution", xlab = "Value",
                              savename = NULL, height = 5, width = 7) {

  if (check_normal && (length(unique(values)) != 1)) {
    Rubrary::check_normal(values)
  }

  if(ctr_measure){
    mea <- mean(as.numeric(values), na.rm = T)
    med <- stats::median(as.numeric(values), na.rm = T)
    sbtt <- paste0("Median: ", round(med, 2), "; Mean: ", round(mea, 2))
  } else {
    sbtt <- ggplot2::waiver()
  }

  plt <- ggplot(data.frame(x = values), aes(x = x)) +
    {if(hist) geom_histogram(aes(y = after_stat(density)), color = "black", fill = "white")} +
    geom_density(alpha = 0.2, fill = "red") +
    {if(ctr_measure) geom_vline(aes(xintercept = mea), color = "blue", linetype = "dashed")} +
    {if(ctr_measure) geom_vline(aes(xintercept = med), color = "red", linetype = "dashed")} +
    labs(
      title = title,
      subtitle = sbtt,
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

#' Density plot by group
#'
#' @import ggplot2
#'
#' @param df dataframe
#' @param value string; colname of values
#' @param group string; colname of group info
#' @param group2 string; group of interest
#' @param pval logical; get KSpval first group vs others
#' @param colors vector; # of colors = to # of groups
#' @param savename string; filepath to save PNG under
#' @param title string; plot title
#' @param xlab string; x-axis label
#' @param rug logical; put rug plot on x axis
#' @param rug_label logical; label rug ticks with first df column
#' @param rug_lab_nudge numeric; fine-tune distance of label from x-axis
#' @param width numeric; ggsave width of plot
#' @param height numeric; ggsave height of plot
#'
#' @return Density plot by group
#' @export
#'
#' @examples
#' set.seed(13)
#' df = data.frame(
#'   group = c(rep("A", 50), rep("B",50)),
#'   values = c(rnorm(50, mean = 0), rnorm(50, mean = 2)))
#' plot_density(df, value = "values", group = "group")

plot_density <- function(df, value, group, group2 = NA, title = NA, pval = T,
                         xlab = value, colors = c("firebrick3", "gray"),
                         rug = F, rug_label = F, rug_lab_nudge = 0.00001,
                         savename = NA, width = 6, height = 4) {
  if(is.na(group2)) {
    group2 <- unique(df[,group])[1]
  }

  # KS pval
  if (pval) {
    ks_pval <- Rubrary::get_kspval(df, value, group, group2)
  }

  dplt <- ggplot(data = df,
                 aes(x = .data[[value]], fill = .data[[group]])) +
    {if(rug) geom_rug(data = df[df[,group] == group2,], color = colors[1])} +
    geom_density(alpha = 0.85) +
    scale_fill_manual(values = colors) +
    {if(rug_label) ggrepel::geom_text_repel(
      data = df[df[,group] == group2,],
      mapping = aes(x = .data[[value]], y = 0, label = df[df[,group] == group2,][,1],),
      angle = 90, hjust = 0, direction = 'x', nudge_y = rug_lab_nudge)} +
    ylab("Density") +
    xlab(xlab) +
    {if (!is.na(title)) labs(title = title)} +
    {if (pval) labs(subtitle = paste0("KS enrich. p-value = ", signif(ks_pval, digits = 4)))} +
    theme_classic()

  if (!is.na(savename)) {
    ggsave(
      filename = savename,
      plot = dplt,
      width = width, height = height
    )
  }
  return(dplt)
}

#' Plot simple scatter
#'
#' Correlation value and method gets placed in plot caption / bottom-right.
#'
#' @import ggplot2
#'
#' @param df dataframe; includes both signatures
#' @param xval string/numeric vector; colname for x axis
#' @param yval string/numeric vector; colname for y axis
#' @param label string; colname for point labels
#' @param group string; colname for group color
#' @param colors char vector; list of colors, where length = # unique groups
#' @param title string; plot title
#' @param subtitle string; plot subtitle
#' @param rank logical; change metric to rank
#' @param xlabel string; xval description
#' @param ylabel string; yval description
#' @param cormethod string; correlation method for stats::cor() or "none"
#' @param pt_size numeric; point size
#' @param pt_alpha numeric; point alpha value
#' @param lbl_size numeric; label text size
#' @param density logical; show density plot along both axes
#' @param guides logical; T to include fitted linear model line
#' @param reverse logical; T to reverse both axes
#' @param savename string; filepath to save figure under
#' @param width numeric; plot width
#' @param height numeric; plot height
#'
#' @return Simple scatter plot as ggplot2
#' @export
#'
#' @examples
#' df <- data.frame(A = 1:10, B = 11:20, C = LETTERS[1:10], D = rep(LETTERS[1:5], each = 2))
#' plot_scatter(df, xval = "A", yval = "B", label = "C", group = "D", title = "Example Scatter")
#'
plot_scatter <- function(df = NULL, xval, yval, label = NULL, group = NULL, colors = NULL, rank = FALSE,
                         xlabel = ifelse(methods::is(xval, "character"), xval, "X"),
                         ylabel = ifelse(methods::is(yval, "character"), yval, "Y"),
                         cormethod = c("pearson", "spearman", "none"), guides = TRUE,
                         pt_size = 2, pt_alpha = 1, lbl_size = 3, density = FALSE,
                         reverse = FALSE, title = NULL, subtitle = NULL,
                         savename = NULL, width = 8, height = 8) {
  cormethod <- match.arg(cormethod)

  if(is.null(df)){ # Create dataframe from vectors if not passed as df
    if (!is.null(label)){
      df <- data.frame(xval = xval, yval = yval, label = label)
    } else {
      df <- data.frame(xval = xval, yval = yval)
    }
    xval <- "xval"
    yval <- "yval"
    # label <- "label"
  } else {
    df <- df[, c(xval, yval, label, group)]
  }

  if(is.null(group)){ # No group
    group <- "group"
    df$group <- "none"
    if(is.null(colors)){ colors = "black" }
  } else { # Group exists
    if(is.null(colors)){ colors = scales::hue_pal()(length(unique(df[[group]])))}
  }

  if (rank) {
    df <- df[order(df[,xval], decreasing = T),]
    df[,xval] <- 1:nrow(df)
    df <- df[order(df[,yval], decreasing = T), ]
    df[,yval] <- 1:nrow(df)
    cormethod <- "spearman"
    # corcoef <- "rho"
    corcoef <- "\U03C1"
    guides <- FALSE
    xlabel <- paste0(xlabel, " Rank")
    ylabel <- paste0(ylabel, " Rank")
  } else {
    corcoef <- "R"
  }

  if(cormethod != "none"){
    corr <- signif(cor(
      x = df[, xval],
      y = df[, yval],
      method = cormethod
    ), 2)
  }

  plt <- ggplot(df, aes(x = .data[[xval]], y = .data[[yval]], color = .data[[group]])) +
    {if(guides) geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red")} +
    geom_point(alpha = pt_alpha, size = pt_size) +
    scale_color_manual(values = colors) +
    xlab(xlabel) +
    ylab(ylabel) +
    {if(reverse) scale_x_reverse()} +
    {if(reverse) scale_y_reverse()} +
    labs(title = title, subtitle = subtitle) +
    {if(cormethod != "none") labs(caption = paste0(
      corcoef, " = ", corr, "; ", tools::toTitleCase(cormethod), " correlation"))} +
    {if(!is.null(label)) ggrepel::geom_text_repel(
      aes(label = .data[[label]]), color = "black", max.overlaps = Inf, size = lbl_size)} +
    theme_classic() +
    {if(length(unique(df[[group]])) == 1) theme(legend.position = "none")}

  if (!is.null(savename)) {
    ggplot2::ggsave(
      filename = savename,
      plot = plt,
      height = height,
      width = width
    )
  }

  if(density) {
    Rubrary::use_pkg("ggExtra")
    mplt <- ggExtra::ggMarginal(
      plt, y = "Density", type = "density", margins = "both",
      size = 6, groupColour = !is.null(group), groupFill = !is.null(group)
    )
    if (!is.null(savename)) {
      savename <- paste0(tools::file_path_sans_ext(savename), "_density.", tools::file_ext(savename))
      ggsave(
        plot = mplt,
        filename = savename,
        height = height, width = width
      )
    }
    plt <- mplt
  }
  return(plt)
}
