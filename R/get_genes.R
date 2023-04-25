utils::globalVariables(c(
  "hgnc_symbol", "gene"
))

#' Retrieve protein coding genes from Ensembl BioMart
#'
#' @return character vector of protein coding genes
#' @export
#'
#' @examples
#' head(Rubrary::get_PC_genes())
get_PC_genes <- function(){
  Rubrary::use_pkg("biomaRt")
  # Select Homo sapiens mart
  mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",
                              mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL"))
  PC_genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "gene_biotype"),
                             filters = "biotype",
                             values = "protein_coding",
                             mart = mart) %>%
    dplyr::filter(hgnc_symbol != "") %>%
    dplyr::rename(gene = hgnc_symbol) %>%
    dplyr::select(gene)

  return(PC_genes$gene)
}

#' Filter dataframe to protein coding genes
#'
#' @param df dataframe; rownames must be gene symbols
#' @param PC_genes character vector; protein coding genes
#' @param search logical; T to use `Seurat::UpdateSymbolList` to match gene symbols better
#'
#' @return dataframe with gene rownames filtered to protein coding only
#' @export
filter_PC_genes <- function(df, PC_genes = Rubrary::get_PC_genes(),
                            search = FALSE){
  nomatch <- rownames(df)[!(rownames(df) %in% PC_genes)]
  warning(paste0("Genes not found in gene list: \n",
                 paste(nomatch, collapse = ", ")),
          call. = FALSE, immediate. = TRUE)

  # Seurat symbol update if desired
  if ((search == FALSE) && (length(nomatch) > 0)){
    if (utils::menu(c("Yes", "No"),
                    title = "\nSearch gene symbols?") == "1") {
      search = TRUE
    }
  }
  if(search){
    Rubrary::use_pkg("Seurat", strict = TRUE)
    updated <- Seurat::UpdateSymbolList(symbols = nomatch)
    names(updated) <- nomatch
    for(g in names(updated)){
      PC_genes[PC_genes == updated[g]] <- g
    }
  }

  # Filter gene rownames
  df_pc <- df[rownames(df) %in% PC_genes,,drop = F]
  return(df_pc)
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
#' get_gene_desc(genes)
get_gene_desc <- function(genes,
                          attrs = c("hgnc_symbol", "description"),
                          clean_desc = TRUE, mart = NULL, verbose = FALSE){
  Rubrary::use_pkg("biomaRt")
  # Select Homo sapiens mart
  if(is.null(mart)){
    mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",
                                mart = biomaRt::useMart(
                                  "ENSEMBL_MART_ENSEMBL", verbose = verbose),
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
