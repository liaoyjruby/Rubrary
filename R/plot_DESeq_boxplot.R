#' Ashvath's DESeq boxplot code in function format
#'
#' @param sig DESeq signature
#' @param data log2UQ transformed data
#' @param anno annotation file
#' @param batch short descriptor for set that will get appended to filename
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter scale_color_manual theme_classic labs
#' @importFrom rlang .data
plot_DESeq_boxplot <- function(sig, data, anno, batch){
  require(ggplot2)
  require(gridExtra)
  require(tibble)

  if(.Platform$OS.type == "unix"){dir.db = "/Users/liaoyj/Dropbox"} else {dir.db = "C:/Users/rubsy/Dropbox"}
  signature <- sig
  merged_data <- data
  merged_anno <- anno
  merged_anno <- tibble::rownames_to_column(merged_anno, "Sample")
  # Get top/bot n=5 genes
  pos_genes <- signature[1:5,1]
  neg_genes <- signature[(nrow(signature)-4):nrow(signature),1]
  top_genes <- pos_genes
  bottom_genes <- neg_genes
  pos_genes <- merged_data[merged_data$gene %in% pos_genes,]
  neg_genes <- merged_data[merged_data$gene %in% neg_genes,]
  pos_genes <- pos_genes[match(top_genes, pos_genes$gene),]
  neg_genes <- neg_genes[match(bottom_genes, neg_genes$gene),]
  # Set row names to genes and rm gene column
  row.names(pos_genes) <- pos_genes$gene
  pos_genes <- pos_genes[-1]
  row.names(neg_genes) <- neg_genes$gene
  neg_genes <- neg_genes[-1]
  # Transpose
  pos_genes <- as.data.frame(t(pos_genes))
  pos_genes <- data.frame(Sample = rownames(pos_genes), pos_genes)
  rownames(pos_genes) <- c() # Rm rownames
  neg_genes <- as.data.frame(t(neg_genes))
  neg_genes <- data.frame(Sample = rownames(neg_genes), neg_genes)
  rownames(neg_genes) <- c()

  # Merge with annotation file to get Batch (UCLA/Patch) and Condition (Primary/Recurrent) info
  pos_input <- merge(pos_genes, merged_anno)
  neg_input <- merge(neg_genes, merged_anno)
  pos_input <- pos_input[c(1:8)]
  neg_input <- neg_input[c(1:8)]
  colnames <- c("Sample", "gene1", "gene2", "gene3", "gene4", "gene5","Batch", "Type")
  pos_gene_names <- colnames(pos_input)[2:6]
  neg_gene_names <- colnames(neg_input)[2:6]
  colnames(pos_input) <- colnames
  colnames(neg_input) <- colnames

  # Plot positive
  pos_gene_1 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene1)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Pos ", pos_gene_names[1]),
         y = "Gene Exp")
  pos_gene_2 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene2)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Pos ", pos_gene_names[2]),
         y = "Gene Exp")
  pos_gene_3 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene3)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Pos ", pos_gene_names[3]),
         y = "Gene Exp")
  pos_gene_4 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene4)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Pos ", pos_gene_names[4]),
         y = "Gene Exp")
  pos_gene_5 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene5)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Pos ", pos_gene_names[5]),
         y = "Gene Exp")

  # Plot negative
  neg_gene_1 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene1)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Neg ", neg_gene_names[1]),
         y = "Gene Exp")
  neg_gene_2 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene2)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Neg ", neg_gene_names[2]),
         y = "Gene Exp")
  neg_gene_3 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene3)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Neg ", neg_gene_names[3]),
         y = "Gene Exp")
  neg_gene_4 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene4)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Neg ", neg_gene_names[4]),
         y = "Gene Exp")
  neg_gene_5 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene5)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col=Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = c("goldenrod2", "navyblue")) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = paste0("Neg ", neg_gene_names[5]),
         y = "Gene Exp")

  # Combine to png
  posGrid <- gridExtra::arrangeGrob(
    pos_gene_1,
    pos_gene_2,
    pos_gene_3,
    pos_gene_4,
    pos_gene_5
  )
  ggplot2::ggsave(filename = paste0(dir.db, "/OV/DESeq/OV_UCLA&Patch_recur_DESeq_log2UQ_Pos_", batch,".png"),
         plot = posGrid, units = "in", width = 13, height = 9, dpi = 300)

  negGrid <- gridExtra::arrangeGrob(
    neg_gene_1,
    neg_gene_2,
    neg_gene_3,
    neg_gene_4,
    neg_gene_5
  )
  ggplot2::ggsave(filename = paste0(dir.db, "/OV/DESeq/OV_UCLA&Patch_recur_DESeq_log2UQ_Neg_", batch,".png"),
         plot = negGrid, units = "in", width = 13, height = 9, dpi = 300)

}
