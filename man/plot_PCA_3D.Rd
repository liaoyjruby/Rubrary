% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCA.R
\name{plot_PCA_3D}
\alias{plot_PCA_3D}
\title{Plot 3D PCA scores/loadings via \code{plotly}}
\usage{
plot_PCA_3D(
  df_pca,
  PCs = c(1:3),
  type = c("Scores", "Loadings"),
  anno = NULL,
  annoname = "Sample",
  annotype = "Type",
  annolabel = annoname,
  label = FALSE,
  colors = NULL,
  title = NULL,
  savename = NULL,
  rotate = FALSE,
  width = 10,
  height = 10
)
}
\arguments{
\item{df_pca}{string or \code{prcomp} obj; (path to) PCA output}

\item{PCs}{num vector; list of numeric PCs to plot (ex. \code{c(1:3)} for first 3 PCs)}

\item{type}{\code{c("Score", "Loading")}}

\item{anno}{string/df; Annotation info for \code{df_pca} with \code{annoname}, \code{annotype}, and \code{annolabel} columns}

\item{annoname}{string; Colname in \code{anno} matching point name}

\item{annotype}{string; Colname in \code{anno} with info to color by}

\item{annolabel}{string; Colname in \code{anno} to label points by, defaults to \code{annoname}}

\item{label}{logical; T to label points}

\item{colors}{char vector; For discrete \code{annotype}, length should be number of unique \code{annotype}s or length 2 for continuous \code{annotype}s where \code{colors[1]} represents low values and \code{colors[2]} represents high values.}

\item{title}{string; Plot title}

\item{savename}{string; File path to save plot under, \code{html} if not an image format}

\item{rotate}{logical; T to have the HTML widget automatically rotate when opened. Only applicable if \code{savename} is not \code{NULL}}

\item{width}{numeric; Saved plot width}

\item{height}{numeric; Saved plot height if saving as image format}
}
\value{
\code{plotly} object
}
\description{
Image saving requires \code{kaleido} python package setup via \code{reticulate} R package. See \code{?plotly::save_image} for more details and setup instructions. Pass in \code{savename} with \code{html} file extension to avoid these errors.
}
\examples{
data(iris)
iris$Sample = rownames(iris)
plot_PCA_3D(
  df_pca = Rubrary::run_PCA(t(iris[,c(1:4)]), screeplot = FALSE),
  PCs = c(1:3),
  type = "Scores",
  anno = iris[,c("Sample", "Species")],
  annoname = "Sample", annotype = "Species",
  title = "Iris PCA Scores 3D"
)

}
