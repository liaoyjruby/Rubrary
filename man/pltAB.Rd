% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scRNA_hypergeo_coexp.R
\name{pltAB}
\alias{pltAB}
\title{Plot gene co-expression scatter for genes A vs B}
\usage{
pltAB(geneB, geneA, hyp_df, sobj, group)
}
\arguments{
\item{geneB}{string; gene B}

\item{geneA}{string; gene A}

\item{hyp_df}{dataframe; hypergeometric p-value table}

\item{sobj}{Seurat object}

\item{group}{string; metadata value to group by}
}
\value{
Scatter plot
}
\description{
Helper function for 'genecoexp_scatter_hyper'
}
\keyword{internal}
