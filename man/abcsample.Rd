% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc.R
\name{abcsample}
\alias{abcsample}
\title{generate ABC samples for causal variant models}
\usage{
abcsample(Zo, snp.data, weights, freq, eps, outdir, N0, N1, nsam_model = 100,
  nsam_gamma = 1000, ncv = 1, uniqify = TRUE, method = c("random",
  "gapfill"))
}
\arguments{
\item{Zo}{observed Z score}

\item{snp.data}{SnpMatrix of reference data}

\item{weights}{LDAK weights}

\item{freq}{phased reference haplotypes with Probability column}

\item{eps}{maximum difference between simulated and observed data
for samples to pass}

\item{outdir}{directory within which output file will be created}

\item{N0}{number of controls}

\item{N1}{number of cases}

\item{nsam_model}{number of causal models to consider}

\item{nsam_gamma}{number of gamma to consider per model}

\item{ncv}{number of causal variants}

\item{uniqify}{if TRUE, don't sample the same model more than nsam_gamma times}

\item{method}{select models at random, or to fill gaps in already sampled}
}
\value{
filename within which samples that passed the threshold
    are written
}
\description{
.. content for \description{} (no empty lines) ..
}
\details{
.. content for \details{} ..
}
\author{
Chris Wallace, Mary Fortune
}
