#' Sum Z Score
#'
#' Slightly modified from Favour's function for OV project
#'
#' @param df numeric dataframe (e.g gene expression df); rownames should contain name of genes/features in signature
#' @param sig vector/list of genes (if applicable, list should be ordered so enriched genes are at the top & de-enriched genes are at the bottom)
#' @param top.and.bottom default = F; if gene list contain both enriched & de-enriched/depleted genes, set to T
#' @param n number of genes used. if top.and.bottom == T, the top n & bottom n genes would be used
#'
#' @return scores samples using signature genes
#' @export
#'
calc_sumZscore <- function(df, sig, topn = length(sig), botn = NA){

    top.genes = sig[1:topn]
    df.top = df[rownames(df) %in% top.genes,]
    print(paste0("Overlapping Top Genes = ", nrow(df.top), " / ", topn))
    df.top.sumz = as.data.frame(
      cbind(
        colnames(df.top),
        rowSums(
          scale(
            t(df.top),
            center = T,
            scale = T
          ),
          na.rm = T
        )
      ),
      stringsAsFactors = F
    )
    colnames(df.top.sumz) = c("sample", "sum.z.score.top")

  if(!is.na(botn)){
    bottom.genes = utils::tail(sig, botn)
    df.bot = df[rownames(df) %in% bottom.genes,]
    print(paste0("Overlapping Bot Genes = ", nrow(df.bot), " / ", botn))
    df.bot.sumz = as.data.frame(
      cbind(
        colnames(df.bot),
        (rowSums(
          scale(
            t(df.bot),
            center = T,
            scale = T
          ),
          na.rm = T
        ) * -1)
      ),
      stringsAsFactors = F
    )
    colnames(df.bot.sumz) = c("sample", "sum.z.score.bottom")

    df.sumz = merge(df.top.sumz, df.bot.sumz, "sample")
    df.sumz$sum.z.score = as.numeric(df.sumz[,2]) + as.numeric(df.sumz[,3])
  }else{
    df.sumz = df.top.sumz
    colnames(df.sumz)[2] = "sum.z.score"
  }
  df.sumz$sum.z.score = as.numeric(df.sumz$sum.z.score)
  df.sumz.ord = df.sumz[order(df.sumz$sum.z.score, decreasing = T), c("sample", "sum.z.score")]

  return(df.sumz.ord)
}
