utils::globalVariables(c(
  "hgnc_symbol", "gene"
))

#' Convert list of genes via BioMart
#'
#' See `biomaRT::listAttributes(mart)` for list of all options for `from_to`.
#'
#' @param genes char vector; genes
#' @param from_to length 2 char vector; 1st value is current format (ex. `ensemb_gene_id`), 2nd value is desired format (ex. `hgnc_symbol`)
#' @param mart biomaRt Mart object
#' @param table logical; T to return conversion table
#'
#' @return char vector with converted gene list to desired format
#' @export
convert_genes <- function(genes, from_to = c("ensembl_gene_id", "hgnc_symbol"),
                          mart = NULL, table = FALSE){
  Rubrary::use_pkg("biomaRt", strict = T)
  if(is.null(mart)){
    mart <- biomaRt::useDataset(
      dataset = "hsapiens_gene_ensembl",
      mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                              host = "https://useast.ensembl.org"))
  }
  conv <- biomaRt::getBM(
    filters = from_to[1],
    attributes = from_to,
    values = genes,
    mart = mart
  )

  if(!table) {conv <- conv$hgnc_symbol}

  return(conv)
}

#' Retrieve protein coding genes from Ensembl BioMart
#'
#' @param mart biomaRt Mart object
#'
#' @return character vector of protein coding genes
#' @export
#'
#' @examples
#' # head(Rubrary::get_PC_genes())
#'
get_PC_genes <- function(mart = NULL){
  Rubrary::use_pkg("biomaRt", strict = T)
  if(is.null(mart)){
    mart <- biomaRt::useDataset(
      dataset = "hsapiens_gene_ensembl",
      mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                              host = "https://useast.ensembl.org"))
  }
  PC_genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "gene_biotype"),
                             filters = "biotype",
                             values = "protein_coding",
                             mart = mart) %>%
    dplyr::filter(hgnc_symbol != "") %>%
    dplyr::rename(gene = hgnc_symbol) %>%
    dplyr::select(gene)

  return(PC_genes$gene)
}

#' Filter dataframe by list of (protein coding) genes
#'
#' @param df dataframe; gene names in `gene_col` or in rownames of dataframe
#' @param gene_col string; colnames of genes, assumed `rownames(df)` if NULL
#' @param genes_filt char vector; (protein coding) genes to include / filter by
#' @param search logical; `TRUE` to use `Seurat::UpdateSymbolList` to match gene symbols better
#'
#' @return dataframe with gene rownames filtered to protein coding only
#' @export
#'
#' @examples
#' df = data.frame(
#'   gene = c("geneA", "geneB", "geneC", "geneD"),
#'   sampA = c(1, 2, 3, 4),
#'   sampB = c(2, 3, 4, 5),
#'   sampC = c(3, 4, 5, 6)
#' )
#' df_filt <- Rubrary::filter_genes(df, genes_filt = c("geneB", "geneC"), gene_col = "gene")
#'
filter_genes <- function(df, genes_filt = Rubrary::get_PC_genes(),
                         gene_col = NULL, search = FALSE){
  if(is.character(df)){
    savename <- tools::file_path_sans_ext(df)
    df <- Rubrary::rread(df, row.names = 1)
  } else {
    savename <- NULL
  }
  if(is.null(gene_col)){
    genes_orig <- rownames(df)
  } else {
    genes_orig <- df[,gene_col]
  }
  nomatch <- genes_filt[!(genes_filt %in% genes_orig)]
  if (length(nomatch) > 0) {
    if(length(nomatch) > 100){
      msg <- paste0(length(nomatch), " genes not found in dataframe gene list: \n",
                    paste(nomatch[1:10], collapse = ", "), ", ...\n",
                    "** Search may take a long time to run!")
    } else {
      msg <- paste0("Genes not found in dataframe gene list: \n",
             paste(nomatch, collapse = ", "), "\n** Total: ",
             length(nomatch))
    }
    warning(msg, call. = FALSE, immediate. = TRUE)
  }

  # Seurat symbol update if desired
  if (search && (length(nomatch) > 0)){
    if (utils::menu(c("Yes", "No"), title = "\nSearch gene symbols?") == "1") {
      search = TRUE
    } else {
      message("** Search gene name match not performed!")
    }
  }
  if(search){
    Rubrary::use_pkg("Seurat", strict = TRUE)
    updated <- Seurat::UpdateSymbolList(symbols = nomatch)
    names(updated) <- nomatch
    for(g in names(updated)){
      genes_filt[genes_filt == updated[g]] <- g
    }
  }

  if(is.null(gene_col)){
    df_filt <- df[rownames(df) %in% genes_filt,,drop = F]
  } else {
    df_filt <- df[df[,gene_col] %in% genes_filt,,drop = F]
  }

  if(!is.null(savename)){
    Rubrary::rwrite(
      x = df_filt, file = paste0(savename, "_filt.txt")
    )
  }

  return(df_filt)
}

#' Get gene description from BioMart
#'
#' @param genes character vector; genes of interest
#' @param attrs character vector; valid attributes from `biomaRt::listAttributes(mart)`
#' @param clean_desc logical; T to remove "\[Source: ...\]" from description field
#' @param verbose logical; output `biomaRt` messages to console
#' @param mart Mart object; will generate one if not provided
#'
#' @return dataframe with attrs columns, populated from BioMart
#' @export
#'
#' @examples
#' genes <- c("POU5F1", "SOX2", "KLF4", "MYC")
#' \donttest{
#'   head(Rubrary::get_gene_desc(genes))
#' }
#'
get_gene_desc <- function(genes,
                          attrs = c("hgnc_symbol", "description"),
                          clean_desc = TRUE, mart = NULL, verbose = FALSE){
  Rubrary::use_pkg("biomaRt")
  # Select Homo sapiens mart
  if(is.null(mart)){
    mart <- biomaRt::useDataset(
      dataset = "hsapiens_gene_ensembl",
      mart = biomaRt::useMart(
        "ENSEMBL_MART_ENSEMBL",
        host="https://useast.ensembl.org",
        verbose = verbose),
      verbose = verbose)
  }

  # Run biomaRt query
  tbl <- biomaRt::getBM(
    attributes = attrs,
    filters  = "hgnc_symbol",
    values = genes,
    mart = mart,
    verbose = verbose)
  if(clean_desc){
    tbl$description <- sub("\\[Source:.*", "", tbl$description)
  }
  return(tbl)
}
