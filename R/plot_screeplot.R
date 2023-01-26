
plot_screeplot <- function(obj_prcomp){
    var_explained = pca$sdev^2 / sum(pca$sdev^2)
    df_scrplt <- data.frame(PC = colnames(pca$rotation),
                            Var.Exp = var_explained[1:length(colnames(pca$rotation))],
                            Cum.Var.Exp = cumsum(var_explained)[1:length(colnames(pca$rotation))])
    
    df_scrplt <- reshape2::melt(df_scrplt, id.var = "PC")
    df_scrplt$value <- df_scrplt$value * 100
    
    scrplt <- ggplot(data = df_scrplt, aes(x = PC, y = value, col = variable, group = variable)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(breaks = seq(0, 100, by = 10)) +
      xlab("Principal Component") +
      ylab("Variance Explained (%)") +
      labs(title = paste0(basename(savename), " Screeplot")) +
      theme_bw() +
      theme(legend.title = element_blank())
    
    plot(scrplt)
    
    ggsave(
      filename = paste0(savename,"_prcomp_screeplot.png"),
      plot = scrplt,
      width = 6, height = 4
    )
}