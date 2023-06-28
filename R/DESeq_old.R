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
#' Not generalizable to non OV project related datasets yet...
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
    Rubrary::rwrite(x = resOrdered, file = paste0(savename, ".txt"))
    rank_file <- resOrdered[c("gene", "sign_log_p")]
    Rubrary::rwrite(rank_file, paste0(savename, ".rnk"))
  }
  return(resOrdered)
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
