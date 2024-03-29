% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEgenes.R
\name{plot_volcano}
\alias{plot_volcano}
\title{Plot (DE genes) volcano plot}
\usage{
plot_volcano(
  df,
  names = "gene",
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 5e-04,
  FCcutoff = 1,
  xlab_high = NULL,
  xlab_low = NULL,
  xlab_size = 7,
  title = NULL,
  subtitle = NULL,
  savename = NULL,
  height = 12,
  width = 16
)
}
\arguments{
\item{df}{df/string; (path to) test statistics, w/ cols \code{names}, \code{x}, \code{y}}

\item{names}{string; colname in \code{df} w/ feature names (genes)}

\item{x}{string; colname for log2 fold changes}

\item{y}{string; colname for nominal/adj p-values}

\item{pCutoff}{numeric; value to set p value cutoff}

\item{FCcutoff}{numeric; value to set log fold change cutoff}

\item{xlab_high}{string; description of high \code{x} values (\verb{\\u2191} for up arrow!)}

\item{xlab_low}{string; description of low \code{x} values (\verb{\\u2193} for down arrow!)}

\item{xlab_size}{numeric; size of x labels}

\item{title}{string; plot title}

\item{subtitle}{string; plot subtitle}

\item{savename}{string; full name of file to save under}

\item{height}{numeric; saved plot height}

\item{width}{numeric; saved plot width}
}
\value{
ggplot2 volcano plot
}
\description{
Wrapper for \code{EnhancedVolcano::EnhancedVolcano()} with additional caption and saving
}
