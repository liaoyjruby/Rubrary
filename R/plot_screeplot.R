utils::globalVariables(c(
  "PC", "value", "variable"
))

#' Plot screeplot from prcomp
#'
#' @import ggplot2
#'
#' @param obj_prcomp prcomp function output object
#' @param savename string; filepath to save PNG under
#' @param label logical; T for % values at points
#' @param npcs integer; # of PCs to include in plot
#' @param cum_var_exp integer; c(1:100), cumulative variance explained threshold
#'
#' @return Screeplot
#' @export
#'
plot_screeplot <- function(obj_prcomp, npcs = ncol(obj_prcomp$x), label = FALSE,
                           cum_var_exp = 80, savename = NULL){
  var_explained = obj_prcomp$sdev^2 / sum(obj_prcomp$sdev^2)
  df_scrplt <- data.frame(PC = colnames(obj_prcomp$rotation),
                          Var.Exp = var_explained[1:length(colnames(obj_prcomp$rotation))],
                          Cum.Var.Exp = cumsum(var_explained)[1:length(colnames(obj_prcomp$rotation))])

  df_scrplt$PC <- sub("PC", "", df_scrplt$PC)
  message(paste0("** Cumulative var. exp. >= ", cum_var_exp,
                 "% at PC ",df_scrplt[df_scrplt$Cum.Var.Exp >= (cum_var_exp/100),]$PC[1],
                 " (", round(df_scrplt[df_scrplt$Cum.Var.Exp >= (cum_var_exp/100),]$Cum.Var.Exp[1] * 100, 1), "%)"))

  df_scrplt <- df_scrplt[1:npcs,]

  df_scrplt$PC <- factor(df_scrplt$PC, levels = unique(df_scrplt$PC))
  df_scrplt <- reshape2::melt(df_scrplt, id.var = "PC")
  df_scrplt$value <- df_scrplt$value * 100
  df_scrplt$label <- paste0(round(df_scrplt$value, 1), "%")

  title_scr <- ifelse(is.null(savename), "PCA", basename(savename))

  scrplt <- ggplot(data = df_scrplt, aes(x = PC, y = value, col = variable, group = variable, label = label)) +
    geom_hline(yintercept = cum_var_exp, color = "red") +
    geom_line() +
    geom_point() +
    scale_y_continuous(breaks = seq(0, 100, by = 10)) +
    xlab("Principal Component") +
    ylab("Variance Explained (%)") +
    labs(title = paste0(title_scr, " Screeplot")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 90)) +
    {if (label) ggrepel::geom_text_repel(size = 3, show.legend = F)}

  if(!is.null(savename)){
    ggsave(
      filename = paste0(savename,"_prcomp_screeplot.png"),
      plot = scrplt,
      width = 8, height = 6
    )
  }

  return(scrplt)
}
