#' Quantile normalization, upper quartile by default
#'
#' In RNA-seq normalization, upper quartile normalization = divide each read count by 75th quantile of read counts in its sample.
#'
#' @param vals numeric df/vector
#' @param q numeric; quartile to normalize by (ex. `q = 75` is 75th quantile / upper quartile)
#' @param target numeric; value to scale by
#' @param min numeric; filter to values above `min` for quantile normalization
#' @param perl logical; `TRUE` to match output of `quartile_norm.pl` commonly used in Graeber Lab RNA-seq analysis
#'
#' @return Quantile normalized numeric vector or dataframe
#' @export
quantile_norm <- function(
    vals, q = 75, target = 1000, min = 1, perl = FALSE){
  if(is.data.frame(vals)){ # Do entire dataframe at once
    vals <- vals %>%
      mutate(across(all_of(names(.)), quantile_norm, perl = perl))
    vals_uq <- vals
  } else {
    vals_filt <- vals[vals >= min] %>% sort()
    if(perl){
      options(digits = 9)
      # Compatibility with Graeber Lab perl quantile_norm.pl - truncates
      quant_perl <- function(vals, quant) {
        len <- length(vals)
        idx <- ((len-1) * quant) + 1 # Perl arrays start idx at 0, R arrays start idx at 1
        qval <- vals[trunc(idx)] # `trunc` replicates Perl's `int` function
        if (idx > trunc(idx)) { # Index not integer value
          qval <- qval + quant * (vals[trunc(idx) + 1] - qval) # Add difference scaled by quant
        }
        return(qval)
      }
      uquant <- quant_perl(vals_filt, q/100)
      vals_uq <- as.numeric(sprintf("%.4f", (vals * (target / uquant))))
    } else {
      uquant <- stats::quantile(vals_filt)[[paste0(q, "%")]]
      vals_uq <- round(vals * (target / uquant), 4)
    }
  }
  return(vals_uq)
}
