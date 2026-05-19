\name{rank_hbic}
\alias{rank_hbic}
\title{HBIC for Penalized Rank Regression}
\usage{
rank_hbic(x, y, beta, n_original, df_thr = 0)
}
\arguments{
\item{x}{Pairwise-difference design matrix.}
\item{y}{Pairwise-difference response.}
\item{beta}{Coefficient vector.}
\item{n_original}{Original sample size before pairwise transformation.}
\item{df_thr}{Threshold for counting nonzero coefficients.}
}
\value{A numeric HBIC value.}
\description{
Computes the HBIC criterion used in the rank-ADMM simulation scripts.
}

