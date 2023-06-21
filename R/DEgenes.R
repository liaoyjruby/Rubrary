utils::globalVariables(c(
  "baseMean", "pvalue", "sign_log_p"
))

#' Plot gene expression boxplots of selected genes (from DESeq signature)
#'
#' Great for sanity checking, to make sure that differentially expressed genes as reported through DE gene analysis metrics (such as DESeq2) are reflected in gene expression data per sample.
#'
#' Providing DE(seq) gene signature (`sig_DE`) is optional; if not provided, `genes` list must be defined. If DE gene info *is* provided, each subplot's title will have the `metric` value appended in parenthesis.
#'
#' @param data df/string; (path to) numeric gene expression matrix, rownames = genes, colnames = samples
#' @param anno df/string; (path to) annotation table, rownames = samples, columns = metadata
#' @param anno_sig string; colname in `anno` to be used as breaks for x-axis
#' @param anno_type string; colname in `anno` to be used to color jitter points by
#' @param sig_DE df/string, optional; (path to) DE(seq) gene results table
#' @param genes char vector; list of genes to plot, must be present in `data`
#' @param metric string; colname in `sig_DE` to rank genes by for auto gene selection, only used if `genes` not provided
#' @param genes_n integer; if `genes` not provided, get top/bottom `genes_n`
#' @param colors char vector; jitter point colors, length == # of unique `anno_type`s
#' @param title string; plot title
#' @param savename string; file path to save plot under
#' @param height numeric; plot height
#' @param width numeric; plot width
#'
#' @return ggplot object, with facet-wrapped boxplot
#' @export
#'
plot_DEgene_boxplot <- function(data, anno, anno_sig, anno_type = NULL, sig_DE, genes = NULL, metric = "sign_log_p", genes_n = 6, colors = NULL, title = NULL, savename = NULL, height = 8, width = 4){
  # DESeq dataframe
  if(!missing(sig_DE)){ # DESeq signature provided
    if(is.character(sig_DE)){ sig_DE <- Rubrary::rread(sig_DE, row.names = 1)}
    de_df <- sig_DE %>%
      select(all_of(metric)) %>%
      arrange(desc(metric)) %>%
      tibble::rownames_to_column("gene")
    # Get genes first (if null)
    if(is.null(genes)){
      de_df <- de_df %>%
        {rbind(utils::head(., genes_n), utils::tail(., genes_n))}
      genes <- de_df$gene
    } else {
      na_genes <- data.frame(
        gene = setdiff(genes, de_df$gene),
        metric = rep(NA, length(setdiff(genes, de_df$gene)))
      ) %>%
        rename_at(vars(colnames(.)), ~colnames(de_df))

      if(nrow(na_genes) != 0){
        message(paste0("Genes not in DE signature: ", paste0(na_genes[,1], collapse = ", ")))
      }

      de_df <- de_df %>%
        filter(gene %in% genes) %>%
        bind_rows(na_genes)
    }
    de_df <- de_df %>%
      mutate(sign = sign(!!sym(metric)))
  } else { # DESeq signature NOT provided
    if(is.null(genes)){ stop("Must provide `genes` if no DE genes signature provided!")}
  }

  # Gene expression & annotation dataframes
  if(is.character(data)){ exp <- Rubrary::rread(data, row.names = 1) } else { exp <- data }
  if(is.character(anno)){ anno <- Rubrary::rread(anno, row.names = 1)}
  # Get expression values for selected genes
  exp_df <- exp %>%
    tibble::rownames_to_column("gene") %>%
    filter(gene %in% genes) %>%
    arrange(match(gene, genes)) # All genes must be in exp_df

  # Simplify annotation
  anno_df <- anno %>%
    tibble::rownames_to_column("sample") %>%
    select(1, all_of(c(anno_sig, anno_type)))
  if(is.null(anno_type)){ # One color
    anno_type <- "type"
    anno_df$type <- anno_type
    if(is.null(colors)){
      colors = "black"
    } else {
      colors = colors[1]
    }
  }

  exp_anno <- exp_df %>%
    tidyr::pivot_longer(!gene, names_to = "sample", values_to = "exp") %>%
    left_join(., anno_df, by = "sample") %>% # Merge anno data
    filter(!is.na(!!sym(anno_sig))) %>% # No NA for x-axis allowed
    as.data.frame()

  if(!missing(sig_DE)){ # Merge in DE info/table
    exp_anno <- exp_anno %>%
      left_join(., de_df, by = "gene") %>% # Merge DESeq data (conditional?)
      mutate(gene = factor(gene, levels = genes)) # Control order of plotting

    # Make facet labels
    gene_labs <- c(paste0(de_df[!is.na(de_df[,metric]), "gene"], " (", round(de_df[!is.na(de_df[,metric]), metric], 2), ")"),
                   de_df[is.na(de_df[,metric]), "gene"])
    names(gene_labs) <- genes
  } else {
    gene_labs <- stats::setNames(nm = genes, genes)
  }

  # Manage colors
  if(is.null(colors)){colors <- scales::hue_pal()(length(unique(exp_anno[,anno_type])))}

  plt <- ggplot(exp_anno, aes(x = .data[[anno_sig]], y = exp)) +
    # geom_boxplot(aes(fill = .data[[metric]]), outlier.shape = NA) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = .data[[anno_type]]), width = 0.05) +
    # scale_fill_gradient(low = "blue", high = "red") +
    scale_color_manual(values = colors) +
    facet_wrap(vars(gene), ncol = 3, scales = "free_y", labeller = labeller(gene = gene_labs)) +
    ylab("Gene Expression") +
    labs(title = title) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    {if(length(colors) == 1) theme(legend.position = "none")}

  if(!is.null(savename)){
    ggsave(
      plot = plt,
      filename = savename,
      height = height, width = width
    )
  }
  return(plt)
}


#' Output DESeq results as dataframe with signed log2 p metric
#'
#' Filters out all `sign_log_p == NA` rows.
#'
#' @import dplyr
#'
#' @param DE_results DESeq results object
#' @param savename string; filepath to save table under (no ext.)
#' @param rank logical; T to output .rnk file compatible with Java GSEA app
#'
#' @return Dataframe with DESeq results
#' @export
#'
output_DESeq <- function(DE_results, savename = NULL, rank = FALSE) {
  res_df <- as.data.frame(DE_results) %>%
    filter(baseMean != 0) %>% # Filter out all where mean was 0
    mutate(sign_log_p = sign(stat)*-log2(pvalue)) %>% # Calc signed log p
    filter(!is.na(sign_log_p)) %>%
    arrange(desc(sign_log_p)) %>%
    tibble::rownames_to_column(var = "gene")

  if(!is.null(savename)){
    write.table(x = res_df,
                file = paste0(savename, ".txt"),
                sep = "\t", row.names = F, quote = F)
  }

  if(rank){
    res_rank <- res_df %>%
      select(c(gene, sign_log_p)) %>%
      mutate(sign_log_p_rank = 1:nrow(res_df))
    if(!is.null(savename)){
      write.table(res_rank[,c("gene", "sign_log_p_rank")],
                  file = paste0(savename, ".rnk"),
                  sep = '\t', quote = F, row.names = F)
    } else {
      res_df <- res_rank
    }
  }
  return(res_df)
}
