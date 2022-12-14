% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_PCA.R
\name{plot_PCA}
\alias{plot_PCA}
\title{Plot PCA with annotation}
\usage{
plot_PCA(
  df,
  dfanno = NA,
  PCx = "PC1",
  PCy = "PC2",
  PCtype = "Score",
  label = T,
  annoname = "Sample",
  annotype = "Batch",
  annotype2 = NA,
  highlight = NA,
  title = paste0("PCA ", PCtype, " Plot - ", annotype),
  subtitle = NA,
  marginal = F,
  save = F,
  savename = "PCA_plot.png"
)
}
\arguments{
\item{df}{string; Path to PC output}

\item{dfanno}{string or df; Annotation info for DF}

\item{PCx}{string; Component on x-axis}

\item{PCy}{string; Component on y-axis}

\item{PCtype}{c("Loading", "Score")}

\item{label}{logical; T to label points}

\item{annoname}{string; Colname in dfanno matching point name}

\item{annotype}{string; Colname in dfanno with info to color by}

\item{annotype2}{string; Colname in dfanno with info to change shape by}

\item{highlight}{char vector; Specific points to shape differently/labels}

\item{title}{string; Plot title}

\item{subtitle}{string; Subtitle for plot}

\item{marginal}{logical; Show density plot along both axes}

\item{save}{logical; T to save to png}

\item{savename}{string; Output plot name}
}
\value{
PCA output plotted with annotation
}
\description{
Plot PCA with annotation
}
