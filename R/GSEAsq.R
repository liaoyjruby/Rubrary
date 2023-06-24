utils::globalVariables(c(
  "Freq", "Term", "pval", "ES", "Category"
))

#' Get GSEA Squared terms/keywords
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param df_GSEA df/string; (path to) GSEA results with `pathway` and `NES` columns
#' @param filt_freq num vector; frequency filter for terms, `filt_freq[1]` is minimum frequency, `filt_freq[2]` is maximum frequency
#' @param signlogp_base integer; log base when calculating signed log_base p metric
#' @param rep0 numeric; value to replace `pval == 0` results with, use `rep0 = 2.2e-16` (rounded `.Machine$double.eps`) to be same as original function
#' @param savename string; path to save output under, no extension
#' @param verbose logical; `TRUE` to output terms as KS p-value is being calculated
#' @param plot logical; `TRUE` to output barplot of significant terms
#' @param plot_pval numeric; include terms with `pval <= plot_pval` in barplot
#' @param plot_fmt string; file extension to save plot as
#' @param seed numeric; randomization seed
#'
#' @return df w/ terms, frequency, (KS) pval, signed enrichment score, and signlogp
#' @export
#'
#' @examples
#' library(dplyr)
#' # Load data
#' deseq_stats <- setNames(
#'   airway_deseq_res[,"sign_log_p"],
#'   airway_deseq_res[,"hgnc_symbol"]
#' )
#' pthwys <- GSEA_pathways
#' # Run (f)gsea
#' gsea_results <- fgsea::fgsea(
#'   pathways = pthwys,
#'   stats = deseq_stats,
#'   eps = 0.0,
#'   minSize = 15,
#'   maxSize  = 500) %>%
#'   arrange(NES)
#' # Get terms
#' head(Rubrary::get_GSEAsq_terms(gsea_results, verbose = FALSE))
#'
get_GSEAsq_terms <- function(
    df_GSEA, savename = NULL,
    filt_freq = c(5, 500), signlogp_base = 10, rep0 = .Machine$double.xmin,
    verbose = TRUE, plot = TRUE, plot_pval = 1e-05, plot_fmt = "png", seed = 13){
  if(is.character(df_GSEA)){ df_GSEA <- Rubrary::rread(df_GSEA)}
  if(!("rank" %in% names(df_GSEA))){ # no ranks provided
    df_GSEA <- df_GSEA %>%
      filter(!is.na(NES)) %>%
      arrange(NES) %>% # default descending
      mutate(rank = 1:nrow(.)) # rank 1 == lowest value
  }

  term_freq_df <- df_GSEA$pathway %>%
    tolower() %>%
    strsplit(., "_") %>%
    unlist() %>%
    table() %>%
    as.data.frame() %>%
    rename(Term = ".") %>%
    arrange(desc(Freq)) %>%
    filter(Freq > filt_freq[1],
           Freq < filt_freq[2]) %>%
    mutate(Term = as.character(toupper(Term)))

  if(verbose){ message(paste0("Total terms: ", nrow(term_freq_df)))}
  GSEA_pw_rnk <- df_GSEA %>%
    select(pathway, NES, rank)

  term_ksp <- function(term){ # per term KS test & ES
    set.seed(seed)
    freq <- term_freq_df[term_freq_df$Term == term, 2]
    if(verbose){
      message(paste0(
        row.names(term_freq_df[term_freq_df$Term == term,])," - ", term,
        " (", freq, ")"))
    }
    GSEA_pw_rnk$term <- factor(ifelse(grepl(term, GSEA_pw_rnk$pathway), TRUE, FALSE), levels = c(TRUE, FALSE))
    # Use `rep0 = 2.2e-16` (rounded `.Machine$double.eps`) to be same as original function
    ksp <- Rubrary::get_kspval(
      df = GSEA_pw_rnk, value = "rank", group = "term", goi = TRUE, rep0 = rep0, signed = TRUE, viz = FALSE)
    # ES: see orig. GSEA squared function that uses https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R
    return(data.frame(Term = term, Freq = freq, pval = ksp$pval, ES = ksp$ES))
  }

  term_ksp_df <- lapply(term_freq_df$Term, term_ksp) %>% # List w/ each element being a term
    bind_rows() %>% # Turn list of rows into df
    arrange(pval) %>% # Sort descending pval
    mutate(signedlogp = sign(ES) * abs(log(pval, signlogp_base)))

  if(plot){
    term_ksp_sig <- term_ksp_df %>%
      filter(pval <= plot_pval) %>%
      arrange(signedlogp)

    term_barp <- ggplot(term_ksp_sig, aes(x = Term, y = signedlogp, fill = Freq)) +
      geom_bar(stat = "identity") +
      scale_x_discrete("Term", limits = term_ksp_sig$Term) +
      labs(title = "GSEA Terms",
           subtitle = paste0("pval cutoff = ", plot_pval)) +
      coord_flip() +
      theme_classic()

    print(term_barp)
  }

  if(!is.null(savename)){
    Rubrary::rwrite(term_ksp_df, paste0(savename, "_GSEAsq_terms_kspvals.txt"))
    if(plot){
      ggsave(
        plot = term_barp,
        filename = paste0(savename, "_GSEAsq_terms_barplot.", plot_fmt),
        width = 8, height = 16
      )
    }
  }
  return(term_ksp_df)
}

#' Run GSEA Squared on GSEA results
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param df_GSEA df/string; (path to) GSEA results with `pathway` and `NES` columns
#' @param categories char vector; list of category names
#' @param cat_terms char vector; list of category terms with each element being keywords separated by "`|`", ex. "CELL_CYCLE|MITOTIC|DNA_REPLICATION"
#' @param rep0 numeric; value to replace `pval == 0` results with, use `rep0 = 2.2e-16` (rounded `.Machine$double.eps`) to be same as original function
#' @param signlogp_base integer; log base when calculating signed log_base p metric
#' @param get_terms logical; `TRUE` to run `Rubrary::get_GSEAsq_terms` to output filtered table of terms and associated statistics
#' @param terms_filt_freq num vector length 2; `get_GSEAsq_terms` arg, `filt_freq[1]` is minimum frequency, `filt_freq[2]` is maximum frequency
#' @param plot_type `c("jitter", "density")`; GSEA squared output plot type
#' @param plot_pval logical; `TRUE` to append p-value to category name. If `ggtext` package is installed, can do fancy markdown formatting
#' @param plot_fmt string; file extension to save plot as
#' @param title string; plot title
#' @param cat_colors char vector; list of colors corresponding to category
#' @param savename string; path to save outputs under (no extension)
#' @param height numeric; output plot height
#' @param width numeric; output plot width
#' @param seed integer; randomization seed
#'
#' @return `GSEAsq` object with `ggplot` object `plot`, df of original GSEA input with categories `pathways`, df of categories KS statistics `categories`, and df of terms KS statistics `terms` if applicable
#' @export
#'
#' @examples
#' library(dplyr)
#' # Load data
#' airway_deseq_res <- Rubrary::airway_deseq_res
#' deseq_stats <- setNames(
#'   airway_deseq_res[,"sign_log_p"],
#'   airway_deseq_res[,"hgnc_symbol"]
#' )
#' pthwys <- Rubrary::GSEA_pathways
#' # Run (f)GSEA
#' gsea_results <- fgsea::fgsea(
#'   pathways = pthwys,
#'   stats = deseq_stats,
#'   eps = 0.0,
#'   minSize = 15,
#'   maxSize  = 500) %>%
#'   arrange(NES)
#' # Run GSEA Squared
#' GSEAsq_terms <- c("METABOLIC","DNA")
#' # Run GSEA squared
#' GSEAsq <- Rubrary::run_GSEA_squared(
#'   df_GSEA = gsea_results,
#'   get_terms = TRUE, verbose = FALSE,
#'   categories = GSEAsq_terms,
#'   cat_terms = GSEAsq_terms,
#'   plot_pval = TRUE,
#'   plot_type = "jitter"
#' )
#' names(GSEAsq) # Various outputs as list
#' GSEAsq$plot

run_GSEA_squared <- function(
    df_GSEA, categories, cat_terms, rep0 = 2.2e-16, signlogp_base = 10,
    get_terms = FALSE, terms_filt_freq = c(5, 500),
    plot_type = c("jitter", "density"), plot_pval = TRUE, title = NULL, cat_colors = NULL,
    savename = NULL, plot_fmt = "png", height = 8, width = 8, seed = 13, verbose = TRUE){

  plot_type = match.arg(plot_type)
  if(is.null(cat_colors)){ cat_colors <- scales::hue_pal()(length(categories))}
  # Load data ----
  if(is.character(df_GSEA)){ df_GSEA <- Rubrary::rread(df_GSEA)}
  df_GSEA <- df_GSEA %>%
    filter(!is.na(NES)) %>%
    arrange(NES) %>% # default descending
    mutate(rank = 1:nrow(.)) # rank 1 == lowest value

  # Get terms ----
  # Can be run independently of GSEAsq but orig. GSEAsq has it included
  if(get_terms){
    terms_df <- Rubrary::get_GSEAsq_terms(
      df_GSEA, savename, seed = seed, rep0 = rep0,
      filt_freq = terms_filt_freq, signlogp_base = signlogp_base,
      plot_fmt = plot_fmt, verbose = verbose)
  }
  cat_df <- data.frame(Category = categories, Term = cat_terms)

  # Get category pathways ----
  cats_df <- lapply(
    cat_df$Category,
    function(cat){
      df_GSEA %>%
        filter(grepl(cat_df[cat_df$Category == cat, "Term"], pathway)) %>%
        mutate(Category = cat)
    }) %>%
    bind_rows() %>%
    arrange(rank)
  # Remaining pathways not captured by categories
  others_df <- df_GSEA %>%
    filter(!(pathway %in% unique(cats_df$pathway))) %>%
    mutate(Category = "Other")

  # Merge cat + non-cat pathways ----
  GSEAsq_df <- bind_rows(cats_df, others_df) %>%
    arrange(rank) %>%
    mutate(Category = factor(Category, levels = c(cat_df$Category, "Other")))

  # Category ks-pvalues ----
  cat_ksp <- function(cat){ # per term KS test & ES
    set.seed(seed)
    freq <- nrow(GSEAsq_df[GSEAsq_df$Category == cat,])
    if(verbose){ message(paste0("** ",cat, " (", freq, ")"))}
    # Filter to relevant category - original function doesn't do this? leaves duplicates in background
    cat_other_df <- GSEAsq_df %>% filter(Category == cat | Category == "Other")
    # Use `rep0 = 2.2e-16` (rounded `.Machine$double.eps`) to be same as original function
    ksp <- Rubrary::get_kspval(
      df = cat_other_df, value = "rank", group = "Category", goi = cat, rep0 = rep0, signed = TRUE, viz = F)
    # ES: see orig. GSEA squared function that uses https://github.com/franapoli/signed-ks-test/blob/master/signed-ks-test.R
    return(data.frame(Category = cat, Freq = freq, pval = ksp$pval, ES = ksp$ES))
  }

  if(verbose){ message("Calculated category KS p-values...")}
  cat_ksp_df <- lapply(categories, cat_ksp) %>% # List w/ each element being a term
    bind_rows() %>% # Turn list of rows into df
    mutate(signedlogp = sign(ES) * abs(log(pval, signlogp_base)),
           sign = case_when( # ifelse(sign(ES) != 1, "\u002b", "\u2212"))
             # Sign symbols opposite of actual sign but better for intuitive understanding?
             sign(ES) == 1 ~ "\u2212", # minus sign == towards left side
             sign(ES) == 0 ~ "0",
             sign(ES) == -1 ~ "\u002b" # plus sign == towards right side
           )
    )

  # Plot ----
  cat_text_obj <- element_text(size = 10, angle = 0, hjust = 1, color = "black")
  ## Plot pval calculations ----
  if(plot_pval){
    Rubrary::use_pkg("ggtext")
    if((requireNamespace("ggtext", quietly = TRUE))){ # ggtext available - fancy!!
      cat_text_obj <- ggtext::element_markdown(
        size = 12, angle = 0, hjust = 1, color = "black", lineheight = 1.2)
      cat_text <- paste0(
        "**", cat_ksp_df$Category, "**<br>",
        "<span style = 'font-size:10pt;color:gray40;'>*p = ", signif(cat_ksp_df$pval, digits = 2), "* (", cat_ksp_df$sign, ")</span>")
    } else { # No ggtext
      cat_text <- paste0(
        cat_ksp_df$Category, "\n",
        "p = ", signif(cat_ksp_df$pval, digits = 2), " (", cat_ksp_df$sign, ")")
    }
  } else {
    cat_text <- cat_ksp_df$Category
  }

  cat_lab <- stats::setNames(nm = cat_ksp_df$Category, cat_text)
  ## Theme settings ----
  if(plot_type == "density"){ # density plot type
    theme_type = theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = cat_text_obj,
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none")
  } else { # jitter plot type
    theme_type = theme(
      legend.position = "none",
      axis.text.y = cat_text_obj
    )
  }
  ## ggplot call ----
  GSEAsq_plt <- ggplot(GSEAsq_df[GSEAsq_df$Category != "Other",], aes(x = rank, color = Category)) +
    {if(plot_type == "jitter") geom_point(aes(y = Category), size = 0.75,
                                          position = position_jitter(seed = seed))} +
    {if(plot_type == "density") geom_density(aes(fill = Category), alpha = 0.5)} +
    {if(plot_type == "density") geom_rug(aes(color = Category), alpha = 0.5)} +
    {if(plot_type == "density") facet_wrap(
      vars(Category), ncol = 1, scales = "free_y", strip.position = "left",
      labeller = labeller(Category = cat_lab))} +
    xlab("Rank") +
    labs(title = title) +
    scale_color_manual(values = cat_colors) +
    {if(plot_type == "jitter") scale_y_discrete(name = NULL, labels = cat_lab)} +
    theme_classic() +
    theme_type

  # Save to file ----
  if(!is.null(savename)){
    Rubrary::rwrite(
      x = GSEAsq_df,
      file = paste0(savename, "_GSEAsq_pathways.txt")
    )
    Rubrary::rwrite(
      x = cat_ksp_df,
      file = paste0(savename, "_GSEAsq_category_kspvals.txt")
    )
    if(plot_fmt == "pdf"){ dev = grDevices::cairo_pdf } else { dev = NULL }
    ggsave(
      plot = GSEAsq_plt,
      filename = paste0(savename, "_GSEAsq_category_", plot_type, "plot.", plot_fmt),
      height = height, width = width,
      device = dev
    )
  }

  # Return GSEAsq object ----
  GSEAsq_output <- list(
    plot = GSEAsq_plt,
    pathways = GSEAsq_df,
    categories = cat_ksp_df
  )
  if(get_terms){ GSEAsq_output$terms <- terms_df}
  class(GSEAsq_output) <- "GSEAsq"
  return(GSEAsq_output)
}

#' Plot categorical percentile rank comparison between two GSEA squared signatures
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param GSEAsq_df1 df/string; "gsea_squared-df.txt" df or path for GSEAsq result 1
#' @param GSEAsq_df2 df/string; "gsea_squared-df.txt" df or path for GSEAsq result 2
#' @param name1 string; descriptor for GSEAsq result 1
#' @param name2 string; descriptor for GSEAsq result 2
#' @param plot_pval logical; `TRUE` to include KS p-value per category in plot
#' @param title string; overall plot title
#' @param colors vector; first color to signify result 1, second for result 2
#' @param rug logical; TRUE for rugplot below density plot
#' @param savename string; filepath to save outputs under (no extension)
#' @param plot_fmt string; file extension to save plot as
#' @param height numeric; plot height
#' @param width numeric; plot width
#'
#' @return facet wrapped ggplot object
#' @export
plot_GSEAsq_density <- function(
    GSEAsq_df1, GSEAsq_df2, name1, name2, plot_pval = TRUE,
    title = "GSEA Squared", colors = c("firebrick", "lightslateblue"), rug = TRUE,
    savename = NULL, plot_fmt = "png", height = 8, width = 8){
  # Load & prep data ----
  if(is.character(GSEAsq_df1)){ GSEAsq_df1 <- Rubrary::rread(GSEAsq_df1) }
  if(is.character(GSEAsq_df2)){ GSEAsq_df2 <- Rubrary::rread(GSEAsq_df2) }
  ## Add ranks & category if missing ----
  if(!("rank" %in% names(GSEAsq_df1))){
    if("rnk" %in% names(GSEAsq_df1)){
      GSEAsq_df1$rank <- GSEAsq_df1$rnk
    } else {
      GSEAsq_df1$rank <- 1:nrow(GSEAsq_df1)
    }
  }
  if(!("rank" %in% names(GSEAsq_df2))){
    if("rnk" %in% names(GSEAsq_df2)){
      GSEAsq_df2$rank <- GSEAsq_df2$rnk
    } else {
      GSEAsq_df2$rank <- 1:nrow(GSEAsq_df2)
    }
  }
  if(!("Category" %in% names(GSEAsq_df1))){
    if("type" %in% names(GSEAsq_df1)){
      GSEAsq_df1$Category <- GSEAsq_df1$type
    } else {
      stop("GSEAsq_df1: no `category` or `type` column")
    }
  }
  if(!("Category" %in% names(GSEAsq_df2))){
    if("type" %in% names(GSEAsq_df2)){
      GSEAsq_df2$Category <- GSEAsq_df2$type
    } else {
      stop("GSEAsq_df2: no `category` or `type` column")
    }
  }
  ## Calculate percentile rank ----
  GSEAsq_df1 <- GSEAsq_df1 %>%
    mutate(pct_rank = rank / max(rank),
           sig = name1)
  GSEAsq_df2 <- GSEAsq_df2 %>%
    mutate(pct_rank = rank / max(rank),
           sig = name2)
  # Merge GSEAsq signatures ----
  GSEAsq_df <- rbind(GSEAsq_df1, GSEAsq_df2) %>%
    select(pathway, NES, signedlogp, rank, Category, pct_rank, sig) %>%
    filter(!grepl("other", Category, ignore.case = TRUE)) %>%
    mutate(sig = factor(sig, levels = c(name1, name2)),
           Category = droplevels(Category))

  # Calculate category kspvalue ----
  cats <- levels(GSEAsq_df$Category)
  cat_ksp_df <- lapply(
    cats,
    function(cat){
      ksp <- Rubrary::get_kspval(
        df = GSEAsq_df[GSEAsq_df$Category == cat,],
        value = "pct_rank",
        group = "sig",
        goi = name1,
        rep0 = .Machine$double.xmin,
        signed = TRUE
      )
      return(data.frame(Category = cat, pval = ksp$pval, ES = ksp$ES))
    }
  ) %>%
    bind_rows() %>%
    mutate(sign = case_when( # ifelse(sign(ES) != 1, "\u002b", "\u2212"))
             # Sign symbols opposite of actual sign but better for intuitive understanding?
             sign(ES) == 1 ~ "\u2212", # minus sign == towards left side
             sign(ES) == 0 ~ "0",
             sign(ES) == -1 ~ "\u002b" # plus sign == towards right side
           )
    )


  # Plot ----
  cat_text_obj <- element_text(size = 10, angle = 0, hjust = 1, color = "black")
  ## Plot pval calculations ----
  if(plot_pval){
    Rubrary::use_pkg("ggtext")
    if((requireNamespace("ggtext", quietly = TRUE))){ # ggtext available - fancy!!
      cat_text_obj <- ggtext::element_markdown(
        size = 12, angle = 0, hjust = 1, color = "black", lineheight = 1.2)
      cat_text <- paste0(
        "**", cat_ksp_df$Category, "**<br>",
        "<span style = 'font-size:10pt;color:gray40;'>*p = ", signif(cat_ksp_df$pval, digits = 2), "* (", cat_ksp_df$sign, ")</span>")
    } else { # No ggtext
      cat_text <- paste0(
        cat_ksp_df$Category, "\n",
        "p = ", signif(cat_ksp_df$pval, digits = 2), " (", cat_ksp_df$sign, ")")
    }
  } else {
    cat_text <- cat_ksp_df$Category
  }
  cat_lab <- stats::setNames(nm = cat_ksp_df$Category, cat_text)

  ## ggplot call ----
  plt <- ggplot(GSEAsq_df, aes(x = pct_rank)) +
    geom_density(aes(fill = sig, color = sig), alpha = 0.5) +
    {if(rug) geom_rug(aes(color = sig))} +
    scale_fill_manual(values = colors, name = "Signature") +
    scale_color_manual(values = colors, name = "Signature") +
    xlab("Rank (percentile)") +
    ylab(NULL) +
    labs(title = title) +
    theme_classic() +
    facet_wrap(
      vars(Category), ncol = 1, scales = "free_y", strip.position = "left",
      labeller = labeller(Category = cat_lab)) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = cat_text_obj,
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom")
  # Save ----
  if(!is.null(savename)){
    Rubrary::rwrite(
      x = cats_ksp_df,
      paste0(savename, "_GSEAsq_compare_category_kspvals.txt")
    )
    if(plot_fmt == "pdf"){ dev = grDevices::cairo_pdf } else { dev = NULL }
    ggsave(
      plot = plt,
      filename = paste0(savename, "_GSEAsq_compare_densityplot.", plot_fmt),
      height = height, width = width,
      device = dev
    )
  }
  return(plt)
}
