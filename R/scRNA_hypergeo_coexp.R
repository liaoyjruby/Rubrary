#' Hypergeometric test, 0 vs. non-0, on columns of a dataframe
#'
#' Originally used for scRNA gene co-expression analysis
#'
#' @param df dataframe; includes columns "A" and "B"
#' @param A string; colname for values A
#' @param B string; colname for values B
#'
#' @return single dataframe row with A&B names, phyper inputs, pvalues
#' @export
phyper_df <- function(df, A, B){
  df_sub <- df[,c(A, B)]

  x <- nrow(df_sub[(df_sub[,A] != 0) & (df_sub[,B] != 0),])
  m <- nrow(df_sub[(df_sub[,A] != 0),])
  n <- nrow(df_sub[(df_sub[,A] == 0),])
  k <- nrow(df_sub[(df_sub[,B] != 0),])

  # Overenrichment
  pval_overenrich <- phyper(
    q = x-1, # -1 for >= instead of >
    m = m,
    n = n,
    k = k,
    lower.tail = FALSE
  )

  # Underenrichment
  pval_underenrich <- phyper(
    q = x,
    m = m,
    n = n,
    k = k,
    lower.tail = TRUE
  )

  hypervals <- data.frame(
    A = A,
    B = B,
    x_inAinB = x,
    m_inA = m,
    n_noA = n,
    k_inB = k,
    pval_overenrich = pval_overenrich,
    pval_overenrich_sig = ifelse(pval_overenrich < 0.05, "*", ""),
    pval_underenrich = pval_underenrich,
    pval_underenrich_sig = ifelse(pval_underenrich < 0.05, "*", "")
  )
  return(hypervals)
}

#' Plot gene co-expression scatter for genes A vs B
#'
#' Helper function for 'genecoexp_scatter_hyper'
#'
#' @param geneB string; gene B
#' @param geneA string; gene A
#' @param hyp_df dataframe; hypergeometric p-value table
#' @param sobj Seurat object
#' @param group string; metadata value to group by
#'
#' @return Scatter plot
#' @export
pltAB <- function(geneB, geneA, hyp_df, sobj, group) {
  if(geneA != geneB){
    hyp_df <- hyp_df[(hyp_df$A == geneA) & (hyp_df$B == geneB),]
    corr <- hyp_df$corr
    hyp_over <- signif(hyp_df$pval_overenrich, digits = 3)
    hyp_under <- signif(hyp_df$pval_underenrich, digits = 3)

    plt <- Seurat::FeatureScatter(
      sobj,
      feature1 = geneA,
      feature2 = geneB,
      group.by = group,
      shuffle = T
    ) +
      ggplot2::labs(
        title = paste0(geneA, " vs. ", geneB),
        subtitle = paste0("Hyp.p over: ", hyp_over, ", under: ", hyp_under, "; Cor: ", corr)
      )
  } else {
    plt <- NA
  }
  return(plt)
}

#' Plot pairwise gene coexpression scatter and calculate hypergeometric enrichment p-value between two genelists
#'
#' @param seuobj Seurat object
#' @param df_dense dataframe; dense gene expression matrix from Seurat (cells in rows, genes in columns)
#' @param genesA char vector; list of genes
#' @param genesB char vector; list of genes
#' @param group string; metadata value to group by
#' @param savepath string; path to choose from with title sans extension
#' @param ncols integer; number of columns in grid
#'
#' @importFrom stats cor phyper
#' @importFrom utils write.csv
#'
#' @return dataframe of hypergeometric test p-values
#' @export
#'
genecoexp_scatter_hyper <- function(seuobj, df_dense = NULL, genesA, genesB, group = NULL,
                                    savepath, ncols = NULL){
  if(is.null(df_dense)){
    df_dense <- as.data.frame(as.matrix(seuobj@assays$RNA@data))
    df_dense <- as.data.frame(t(df_dense[unique(c(genesA, genesB)),]))
  } else {
    df_dense <- as.data.frame(t(df_dense[unique(c(genesA, genesB)),]))
  }

  df_hypervals <- data.frame()
  for(geneA in genesA) {
    for(geneB in genesB) {
      if(geneA != geneB){
        # debugonce(phyper_df)
        row <- phyper_df(
          df = df_dense,
          A = geneA,
          B = geneB
        )
        corr <- round(cor(x = df_dense[,geneA], y = df_dense[,geneB]), digits = 2)
        row$corr <- corr
        df_hypervals <- rbind(df_hypervals, row)
      } else {
        message("** Gene A == Gene B; skipping...")
      }
    }
    plts_A <- lapply(genesB, pltAB,
                     geneA = geneA,
                     hyp_df = df_hypervals,
                     sobj = seuobj,
                     group = group)

    plts_A <- plts_A[!is.na(plts_A)]

    ncols <- ifelse(is.null(ncols), ifelse(length(plts_A) < 3, length(plts_A), 3), ncols)

    grid_A <- patchwork::wrap_plots(
      plts_A,
      ncol = ncols
    ) + patchwork::plot_layout(guides = "collect") & ggplot2::theme(legend.position = "right")

    ggplot2::ggsave(
      filename = paste0(savepath, "_", geneA, ".png"),
      plot = grid_A,
      width = ncols * 4,
      height = ceiling(length(plts_A)/ncols) * 4
    )
  }

  # Save hypergeometric values
  write.csv(
    x = df_hypervals,
    file = paste0(savepath, "_hypgeo_pval.csv"),
    row.names = F
  )
  write.csv(
    x = df_hypervals[, c(1:2, 7:10)],
    file = paste0(savepath, "_hypgeo_pval_simple.csv"),
    row.names = F
  )

  return(df_hypervals)
}
