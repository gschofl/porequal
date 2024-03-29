% Generated by roxygen2 (4.0.1): do not edit by hand
\name{plot_alignment_quality}
\alias{plot_alignment_quality}
\title{Plot alignment quality.}
\usage{
plot_alignment_quality(x, smooth = TRUE, annotate = TRUE, annot.size = 4)
}
\arguments{
\item{x}{A \code{\link{consensus_table}}.}

\item{annotate}{Annotate mismatch positions.}

\item{annot.size}{If \code{annotate = TRUE}, size of the annotation text.}

\item{smooth}{Plot a smoothing function (gam: \code{y ~ s(x)}).}
}
\value{
Plots the probability to recover the reference sequence and summarises it in a
boxplot.
}
\description{
This function graphs the probability to recover the reference sequence from
an multiple alignment of nanopore reads.
}

