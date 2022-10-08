if (.Platform$OS.type == "unix") {
  dir.db <- "/Users/liaoyj/Dropbox"
} else {
  dir.db <- "C:/Users/rubsy/Dropbox"
}

utils::globalVariables(c(
  "Type", "Batch",
  "gene1", "gene2", "gene3", "gene4", "gene5"
))

#' df_to_mtx
#'
#' @param df A dataframe
#' @param merged If dataframe is made from combined datasets
#' @param filter Apply filter for genes with 9+ zeros
#'
#' @return Dataframe restructured as matrix
#' @export
#'
df_to_mtx <- function(df, merged = T, filter = T) {
  if (!merged) {
    filter <- F
  }
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

run_DESeq <- function(mtx_rawcts, mtx_annotation, outname, merged = F, shrink = F) {
  mtx <- mtx_rawcts
  anno <- mtx_annotation

  if (grepl("UCLA&Patch", outname)) {
    set <- "UCLA&Patch"
  } else if (grepl("UCLA", outname)) {
    set <- "UCLA"
  } else {
    set <- "Patch"
  }

  # Paired
  dds_simple <- DESeq2::DESeqDataSetFromMatrix(
    countData = mtx,
    colData = anno,
    design = ~ Subject + Condition
  )

  if (merged) {
    num_batches <- length(unique(anno$Batch))
    subj_n <- c()
    for (i in 1:num_batches) {
      batch_pairs <- sum(anno$Batch == unique(anno$Batch)[i]) / 2
      subj_n <- c(subj_n, rep(1:batch_pairs, each = 2)) # list of 1-n_batch, together
    }
    anno$Subject.n <- subj_n
    anno$Subject.n <- factor(anno$Subject.n)
    # Single DS DESeq should already have Condition / Subject columns

    design <- stats::model.matrix(~ Batch + Batch:Subject.n + Batch:Condition, anno)
    design <- design[, colSums(design != 0) > 0]

    seq <- DESeq2::DESeq(dds_simple, full = design, betaPrior = FALSE)
    contrast_list <- list(c(DESeq2::resultsNames(seq)[length(DESeq2::resultsNames(seq)) - 1], DESeq2::resultsNames(seq)[length(DESeq2::resultsNames(seq))]))
  } else {
    seq <- DESeq2::DESeq(dds_simple)
    contrast_list <- c("Condition", "Recurrent", "Primary")
  }
  DESeq2::resultsNames(seq)
  res <- DESeq2::results(seq, contrast = contrast_list)
  DESeq2::plotMA(res, ylim = c(-2, 2), main = paste0(set, " DESeq"))
  if (shrink) {
    res_lfc <- DESeq2::lfcShrink(seq, res = res, contrast = contrast_list, type = "ashr")
    DESeq2::plotMA(res_lfc, ylim = c(-2, 2), main = paste0(set, " DESeq log2 Shrink"))
    res_lfc_df <- as.data.frame(res_lfc)
  }
  res_df <- as.data.frame(res)
  res_df$sign_log_p <- sign(res_df$stat) * -log2(res_df$pvalue)
  resOrdered <- res_df[order(res_df$sign_log_p, decreasing = T), ]
  resOrdered <- resOrdered[stats::complete.cases(resOrdered), ] # excludes rows with NA
  resOrdered <- tibble::rownames_to_column(resOrdered, "gene")
  utils::write.table(x = resOrdered, file = paste0(outname, ".txt"), sep = "\t", row.names = F, quote = F)
  rank_file <- resOrdered[c("gene", "sign_log_p")]
  utils::write.table(rank_file, paste0(outname, ".rnk"), sep = "\t", quote = F, col.names = F, row.names = F)
}

#' Filter signature to PC coding only by matching with genes in log2UQcounts
#'
#' @param sig_path Path to DESeq signature
#'
#' @return Text and rnk files subsetted to PC coding only genes
#' @export
#'
make_pconly <- function(sig_path) {
  if (grepl(".rnk", sig_path)) {
    rnkonly <- T
  } else {
    rnkonly <- F
  }
  pathsplit <- strsplit(sig_path, "[_|.]")[[1]]
  set <- pathsplit[2] # UCLA or Patch or UCLA/Patch
  slogp <- pathsplit[5] # slogp or slogponly
  if (set == "Patch") {
    batch <- ""
  } else {
    batch <- paste0("_", pathsplit[length(pathsplit) - 1])
  }
  merged_OV_log2UQ <- utils::read.table(paste0(dir.db, "/OV/OV_UCLA&Patch_recur_log2UQcounts.tsv"), header = T) # Filter to genes that are present in the log2UQ table
  if (!rnkonly) {
    sig <- utils::read.delim(paste0(dir.db, sig_path), header = T)
    sig_pconly <- sig[sig$gene %in% merged_OV_log2UQ$gene, ] # Match with log2UQ coding genes
    utils::write.table(sig_pconly, paste0(dir.db, "/OV/DESeq/OV_", set, "_recur_prevpost_slogp_DESeq_pconly", batch, ".txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  }
  # Make rank
  sig_rnk <- utils::read.delim(paste0(dir.db, "/OV/DESeq/OV_", set, "_recur_prevpost_", slogp, "_DESeq", batch, ".rnk"), header = F, col.names = c("Gene", "Score"))
  sig_pconly_rnk <- sig_rnk[sig_rnk$Gene %in% merged_OV_log2UQ$gene, ] # Match with log2UQ coding genes
  # sig_pconly_rnk$Score_PC <- 1:nrow(sig_pconly_rnk) # Rerank after filter to protein coding only
  utils::write.table(sig_pconly_rnk, paste0(dir.db, "/OV/DESeq/OV_", set, "_recur_prevpost_", slogp, "_DESeq_pconly", batch, ".rnk"), sep = "\t", quote = F, col.names = F, row.names = F)
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
#' @param sig Signature
#' @param data Gene expression data
#' @param anno Annotation file for samples
#' @param batch Descriptive word/phrase appended to output path
#' @param merged If dataframe is made from combined datasets
#'
#' @return Paired boxplots comparing pre-post for top/bot 5 DESeq genes, plotting gene expression data
#' @export
#'
plot_DESeq_boxplot <- function(sig, data, anno, batch, merged = T) {
  if (.Platform$OS.type == "unix") {
    dir.db <- "/Users/liaoyj/Dropbox"
  } else {
    dir.db <- "C:/Users/rubsy/Dropbox"
  }
  signature <- sig
  merged_data <- data
  merged_anno <- anno
  if (!merged) {
    merged_anno <- tibble::rownames_to_column(merged_anno, "Sample")
  }
  # Get top/bot n=5 genes
  pos_genes <- signature[1:5, 1]
  neg_genes <- signature[(nrow(signature) - 4):nrow(signature), 1]
  top_genes <- pos_genes
  bottom_genes <- neg_genes
  pos_genes <- merged_data[merged_data$gene %in% pos_genes, ]
  neg_genes <- merged_data[merged_data$gene %in% neg_genes, ]
  pos_genes <- pos_genes[match(top_genes, pos_genes$gene), ]
  neg_genes <- neg_genes[match(bottom_genes, neg_genes$gene), ]
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
  colnames <- c("Sample", "gene1", "gene2", "gene3", "gene4", "gene5", "Batch", "Type")
  pos_gene_names <- colnames(pos_input)[2:6]
  neg_gene_names <- colnames(neg_input)[2:6]
  colnames(pos_input) <- colnames
  colnames(neg_input) <- colnames

  colors <- c("goldenrod2", "navyblue")

  if (!merged) {
    pos_input$Type <- pos_input$Batch
    pos_input$Batch <- "UCLA"
    neg_input$Type <- neg_input$Batch
    neg_input$Batch <- "UCLA"
    colors <- c("navyblue")
  }

  # Plot positive
  pos_gene_1 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene1)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", pos_gene_names[1]),
      y = "Gene Exp"
    )
  pos_gene_2 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene2)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", pos_gene_names[2]),
      y = "Gene Exp"
    )
  pos_gene_3 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene3)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", pos_gene_names[3]),
      y = "Gene Exp"
    )
  pos_gene_4 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene4)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", pos_gene_names[4]),
      y = "Gene Exp"
    )
  pos_gene_5 <- ggplot2::ggplot(pos_input, ggplot2::aes(Type, gene5)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Pos ", pos_gene_names[5]),
      y = "Gene Exp"
    )

  # Plot negative
  neg_gene_1 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene1)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", neg_gene_names[1]),
      y = "Gene Exp"
    )
  neg_gene_2 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene2)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", neg_gene_names[2]),
      y = "Gene Exp"
    )
  neg_gene_3 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene3)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", neg_gene_names[3]),
      y = "Gene Exp"
    )
  neg_gene_4 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene4)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", neg_gene_names[4]),
      y = "Gene Exp"
    )
  neg_gene_5 <- ggplot2::ggplot(neg_input, ggplot2::aes(Type, gene5)) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(ggplot2::aes(col = Batch), width = 0.05, height = 0) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("Neg ", neg_gene_names[5]),
      y = "Gene Exp"
    )

  if (merged) {
    title <- "OV_UCLA&Patch_recur_DESeq_log2UQ"
  } else {
    title <- "OV_UCLA_recur_DESeq_log2UQ"
  }

  # Combine to png
  posGrid <- gridExtra::arrangeGrob(
    pos_gene_1,
    pos_gene_2,
    pos_gene_3,
    pos_gene_4,
    pos_gene_5
  )
  ggplot2::ggsave(
    filename = paste0(dir.db, "/OV/DESeq/", title, "_Pos_", batch, ".png"),
    plot = posGrid, units = "in", width = 13, height = 9, dpi = 300
  )

  negGrid <- gridExtra::arrangeGrob(
    neg_gene_1,
    neg_gene_2,
    neg_gene_3,
    neg_gene_4,
    neg_gene_5
  )
  ggplot2::ggsave(
    filename = paste0(dir.db, "/OV/DESeq/", title, "_Neg_", batch, ".png"),
    plot = negGrid, units = "in", width = 13, height = 9, dpi = 300
  )
}
