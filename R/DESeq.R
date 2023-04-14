utils::globalVariables(c(
  "Type", "Batch",
  "gene1", "gene2", "gene3", "gene4", "gene5"
))

#' df_to_mtx
#'
#' @param df dataframe
#' @param filter logical; T to filter out genes with 9+ zeros
#'
#' @return Dataframe restructured as matrix
#' @export
#'
df_to_mtx <- function(df, filter = T) {
  mtx <- df
  rownames(mtx) <- df$gene
  mtx <- mtx[, -1]
  mtx <- as.matrix(mtx)
  # mtx <- as.data.frame(mtx)
  mtx <- round(mtx) # Round all data in merged_mtx to integers

  if (filter) {
    # Filter out genes (rows) when there are more than nine zeros for count values
    remove <- c()
    for (i in 1:nrow(mtx))
    {
      if (sum(mtx[i, ] == 0) > 9) {
        remove <- c(remove, i)
      }
    }
    # print(remove)
    mtx_filt <- mtx[-remove, ]
    mtx <- mtx_filt
  }
  return(mtx)
}

#' Run DESeq
#'
#' Not sure if generalizable to non OV project related datasets yet...
#'
#' @param mtx_rawcts matrix; numeric matrix of raw counts
#' @param mtx_annotation matrix; samples as rownames, cols "Subject" & "Condition" (& "Batch" if merged)
#' @param savename string; filename to save outputs under (no ext.)
#' @param merged logical; T if dataset has two batches merged
#' @param shrink logical; apply lfcShrink "ashr" to results (NOT IMPLEMENTED)
#'
#' @return DESeq results in dataframe
#' @export
#'
run_DESeq <- function(mtx_rawcts, mtx_annotation, savename = NULL,
                      merged = F, shrink = F) {
  # Paired
  dds_simple <- DESeq2::DESeqDataSetFromMatrix(
    countData = mtx_rawcts,
    colData = mtx_annotation,
    design = ~ Subject + Condition
  )

  if (merged) {
    num_batches <- length(unique(mtx_annotation$Batch))
    subj_n <- c()
    for (i in 1:num_batches) {
      batch_pairs <- sum(mtx_annotation$Batch == unique(mtx_annotation$Batch)[i]) / 2
      subj_n <- c(subj_n, rep(1:batch_pairs, each = 2)) # list of 1-n_batch, together
    }
    mtx_annotation$Subject.n <- subj_n
    mtx_annotation$Subject.n <- factor(mtx_annotation$Subject.n)
    # Single DS DESeq should already have Condition / Subject columns

    design <- stats::model.matrix(~ Batch + Batch:Subject.n + Batch:Condition, mtx_annotation)
    design <- design[, colSums(design != 0) > 0]

    seq <- DESeq2::DESeq(dds_simple, full = design, betaPrior = FALSE)
    contrast_list <- list(c(DESeq2::resultsNames(seq)[length(DESeq2::resultsNames(seq)) - 1],
                            DESeq2::resultsNames(seq)[length(DESeq2::resultsNames(seq))]))
  } else {
    seq <- DESeq2::DESeq(dds_simple)
    contrast_list <- c("Condition", "Recurrent", "Primary")
  }
  DESeq2::resultsNames(seq)
  res <- DESeq2::results(seq, contrast = contrast_list)
  DESeq2::plotMA(res, ylim = c(-2, 2), main = "DESeq")
  if (shrink) {
    res_lfc <- DESeq2::lfcShrink(seq, res = res, contrast = contrast_list, type = "ashr")
    DESeq2::plotMA(res_lfc, ylim = c(-2, 2), main = "DESeq log2 Shrink")
    res_lfc_df <- as.data.frame(res_lfc)
  }
  res_df <- as.data.frame(res)
  res_df$sign_log_p <- sign(res_df$stat) * -log2(res_df$pvalue)
  resOrdered <- res_df[order(res_df$sign_log_p, decreasing = T), ]
  resOrdered <- resOrdered[stats::complete.cases(resOrdered), ] # excludes rows with NA
  resOrdered <- tibble::rownames_to_column(resOrdered, "gene")
  if(!is.null(savename)){
    utils::write.table(x = resOrdered, file = paste0(savename, ".txt"), sep = "\t", row.names = F, quote = F)
    rank_file <- resOrdered[c("gene", "sign_log_p")]
    utils::write.table(rank_file, paste0(savename, ".rnk"), sep = "\t", quote = F, col.names = F, row.names = F)
  }
  return(resOrdered)
}

#' Filter signature to PC coding only by matching with genes in log2UQcounts
#'
#' @param sig_path string; path to DESeq signature
#'
#' @return DESeq signature subsetted to PC only genes
#' @export
#'
make_pconly <- function(sig_path) {
  if (grepl(".rnk", sig_path)) {
    rnkonly <- T
  } else {
    rnkonly <- F
  }

  pc_genes <- read.delim("/Users/liaoyj/Library/CloudStorage/Dropbox/Ovarian Project/log2_coding_expression_datasets/protein-coding_gene.txt")$symbol
  if (!rnkonly) {
    sig <- utils::read.delim(sig_path, header = T)
    sig_pconly <- sig[sig$gene %in% pc_genes, ] # Match with log2UQ coding genes
    utils::write.table(sig_pconly, paste0(tools::file_path_sans_ext(sig_path), "_pconly.txt"),
                       sep = "\t", quote = F, col.names = T, row.names = F)
  }
  # Make rank
  sig_rnk <- utils::read.delim(paste0(tools::file_path_sans_ext(sig_path), ".rnk"), header = F, col.names = c("gene", "score"))
  sig_pconly_rnk <- sig_rnk[sig_rnk$gene %in% pc_genes, ] # Match with log2UQ coding genes
  utils::write.table(sig_pconly_rnk, paste0(tools::file_path_sans_ext(sig_path), "_pconly.rnk"),
                     sep = "\t", quote = F, col.names = F, row.names = F)

  return(sig_pconly)
}

#' plot_DESeq_scatter
#'
#' @param sig1path Path to 1st signature (x-axis)
#' @param sig2path Path to 2nd signature (y-axis)
#' @param batch Descriptive word/phrase appended to output path
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param title Plot title
#' @param desc Plot subtitle
#' @param rank Plot rank scatter instead of slogp/raw values
#'
#' @return Plot comparing 1 sig to another
#' @export
#'
plot_DESeq_scatter <- function(sig1path, sig2path, batch = "",
                               xlab = "Old DESeq signed log p", ylab = "Updated DESeq signed log p",
                               title = "", desc = "", rank = F) {
  pathsplit <- strsplit(sig2path, "[_|.]")[[1]]
  set <- pathsplit[2] # UCLA or Patch or UCLA/Patch
  if (set == "UCLA&Patch") {
    set2 <- "UCLA & Patch"
  } else {
    set2 <- set
  }

  if (title == "") {
    title <- paste0(set2, " DESeq Signature Old vs Updated")
  }

  sig1 <- utils::read.delim(sig1path, header = T)
  # if (set == "Patch"){
  #   batch <- ""
  #   colnames(sig1) <- c("gene", "sign_log_p")
  # } else if (batch == ""){
  #   batch <- paste0("_",pathsplit[length(pathsplit)-1])
  # }

  sig2 <- utils::read.delim(sig2path, header = T)
  if (rank) {
    sig1 <- sig1[order(sig1$sign_log_p, decreasing = T), ]
    sig2 <- sig2[order(sig2$sign_log_p, decreasing = T), ]
    sig1$sign_log_p <- 1:nrow(sig1)
    sig2$sign_log_p <- 1:nrow(sig2)
  }
  sig_merged <- base::merge(sig1[, c("gene", "sign_log_p")], sig2[, c("gene", "sign_log_p")], by = "gene", all = T)
  colnames(sig_merged) <- c("gene", "sign_log_p.1", "sign_log_p.2")
  sig1_range <- c(min(sig_merged$sign_log_p.1, na.rm = T), max(sig_merged$sign_log_p.1, na.rm = T))
  sig2_range <- c(min(sig_merged$sign_log_p.2, na.rm = T), max(sig_merged$sign_log_p.2, na.rm = T))
  limits <- c(min(sig1_range, sig2_range) * 1.1, max(sig1_range, sig2_range) * 1.1)
  plt <- ggplot2::ggplot(sig_merged, ggplot2::aes_string(x = "sign_log_p.1", y = "sign_log_p.2")) +
    ggplot2::geom_point(alpha = 0.2, size = 0.5) +
    # geom_abline(linetype="dashed", ggplot2::aes(intercept=0, slope=1), size = 1) +
    # geom_smooth(method = "lm", se=FALSE, color = 'red') +
    ggpubr::stat_cor(
      method = "pearson",
      label.x = limits[1] + (limits[2] * 0.05),
      label.y = limits[2]
    ) +
    # stat_regline_equation(label.x = -40, label.y = 28) +
    # coord_fixed(xlim = limits, ylim = limits) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title, # include comp % on axis
      subtitle = desc,
      x = xlab,
      y = ylab
    )
  print(plt)
  ggplot2::ggsave(
    filename = paste0(dirname(sig2path), "/OV_", set, "_recur_DESeq_slogp_oldvnew", batch, ".png"),
    plot = plt,
    height = 8,
    width = 8
  )
}

#' plot_DESeq_boxplot
#'
#' @param sig dataframe; sorted signature with genes as first column
#' @param data dataframe; gene expression data w/ "gene" column
#' @param anno dataframe; sample annotation file w/ "Batch" & "Type" cols
#' @param merged logica; T if dataframe is made from combined datasets
#' @param colors vector; list of colors, n = number of batches
#' @param savename string; filepath to save figures under (no ext.)
#' @param width numeric; width of plot saved
#' @param height numeric; height of plot saved
#' @param batch string; if not merged,
#'
#' @return Boxplots of gene expression for top/bot 5 genes in sig
#' @export
#'
plot_DESeq_boxplot <- function(sig, data, anno, merged = T, batch = "UCLA",
                               colors = c("goldenrod2", "navyblue"),
                               savename = NULL, width = 13, height = 9) {

  # Get top 5 genes
  pos_genes <- sig[1:5, 1]
  top_genes <- pos_genes
  pos_genes <- data[data$gene %in% pos_genes, ]
  pos_genes <- pos_genes[match(top_genes, pos_genes$gene), ]
  rownames(pos_genes) <- NULL
  pos_genes <- tibble::column_to_rownames(pos_genes, var = "gene")
  pos_genes <- as.data.frame(t(pos_genes))
  pos_genes <- tibble::rownames_to_column(pos_genes, var = "Sample")

  # Get bot 5 genes
  neg_genes <- sig[(nrow(sig) - 4):nrow(sig), 1]
  bottom_genes <- neg_genes
  neg_genes <- data[data$gene %in% neg_genes, ]
  neg_genes <- neg_genes[match(bottom_genes, neg_genes$gene), ]
  rownames(neg_genes) <- NULL
  neg_genes <- tibble::column_to_rownames(neg_genes, var = "gene")
  neg_genes <- as.data.frame(t(neg_genes))
  neg_genes <- tibble::rownames_to_column(neg_genes, var = "Sample")

  # Merge with annotation file
  pos_input <- dplyr::left_join(pos_genes, tibble::rownames_to_column(anno, var = "Sample"), by = "Sample")
  neg_input <- dplyr::left_join(neg_genes, tibble::rownames_to_column(anno, var = "Sample"), by = "Sample")
  pos_input <- pos_input[,c("Sample", top_genes, "Batch", "Type")]
  neg_input <- neg_input[,c("Sample", bottom_genes, "Batch", "Type")]
  colnames <- c("Sample", "gene1", "gene2", "gene3", "gene4", "gene5", "Batch", "Type")
  colnames(pos_input) <- colnames
  colnames(neg_input) <- colnames

  if (!merged) {
    pos_input$Batch <- batch
    neg_input$Batch <- batch
  }

  # Plot positive
  pos_gene_1 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene1)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", top_genes[1]),
      y = "Gene Exp"
    )
  pos_gene_2 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene2)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", top_genes[2]),
      y = "Gene Exp"
    )
  pos_gene_3 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene3)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", top_genes[3]),
      y = "Gene Exp"
    )
  pos_gene_4 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene4)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", top_genes[4]),
      y = "Gene Exp"
    )
  pos_gene_5 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene5)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", top_genes[5]),
      y = "Gene Exp"
    )

  # Plot negative
  neg_gene_1 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene1)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", bottom_genes[1]),
      y = "Gene Exp"
    )
  neg_gene_2 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene2)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", bottom_genes[2]),
      y = "Gene Exp"
    )
  neg_gene_3 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene3)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", bottom_genes[3]),
      y = "Gene Exp"
    )
  neg_gene_4 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene4)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", bottom_genes[4]),
      y = "Gene Exp"
    )
  neg_gene_5 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene5)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", bottom_genes[5]),
      y = "Gene Exp"
    )

  # Combine to png
  posGrid <- gridExtra::arrangeGrob(
    pos_gene_1,
    pos_gene_2,
    pos_gene_3,
    pos_gene_4,
    pos_gene_5
  )
  ggplot2::ggsave(
    filename = paste0(savename, "_Pos.png"),
    plot = posGrid, units = "in", width = width, height = height, dpi = 300
  )

  negGrid <- gridExtra::arrangeGrob(
    neg_gene_1,
    neg_gene_2,
    neg_gene_3,
    neg_gene_4,
    neg_gene_5
  )
  ggplot2::ggsave(
    filename = paste0(savename, "_Neg.png"),
    plot = negGrid, units = "in", width = width, height = height, dpi = 300
  )
}
