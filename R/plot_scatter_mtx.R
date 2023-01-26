#
# panel.cor <- function(x, y){
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- round(cor(x, y), digits=2)
#   txt <- paste0("R = ", r)
#   cex.cor <- 0.8/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor * r)
# }
#
# upper.panel<-function(x, y){
#   points(x,y, pch = 19)
# }
#
# plot_scatter_mtx <- function(df, correlation = TRUE){
#   # Correlation panel
#
#
#   # png(filename = paste0(subset, "_TTest.png"), width = 7, height = 7, units = "in", res = 250)
#   if(correlation){
#     pairs(df[,2:7],
#           lower.panel = panel.cor,
#           upper.panel = upper.panel)
#   } else {
#     pairs(df[,2:7],
#           lower.panel = NULL,
#           upper.panel = upper.panel)
#   }
#
# }
#
# plot_PCA_scattermtx <- function(){
#
# }
