\name{simulation_utils}
\alias{make_rank_pairs}
\alias{prepare_rank_data}
\alias{simulate_rank_example}
\alias{rank_admm_metrics}
\title{Rank-ADMM Simulation and Evaluation Utilities}
\usage{
make_rank_pairs(x, y)
prepare_rank_data(x, y, standardize = TRUE, rank_transform = TRUE)
simulate_rank_example(n = 60L, p = 200L, rho = 0.5,
  error = c("normal", "t4", "mixture"), beta = NULL)
rank_admm_metrics(beta_hat, beta_true, threshold = 0)
}
\arguments{
\item{x}{Design matrix.}
\item{y}{Response vector.}
\item{standardize}{Whether to center/scale \code{x} and center \code{y}.}
\item{rank_transform}{Whether to construct pairwise differences.}
\item{n, p, rho}{Simulation sample size, dimension, and AR(1) covariance parameter.}
\item{error}{Error distribution.}
\item{beta}{True coefficient vector.}
\item{beta_hat, beta_true}{Estimated and true coefficient vectors.}
\item{threshold}{Threshold used only for support-size, FP, and FN summaries. The L1 error is always computed from raw coefficients.}
}
\value{Utility-specific objects: transformed data, simulated data, or a metrics data frame.}
\description{
Utilities for reproducing the pairwise rank transformation, Example 1/2 data generation, and raw-beta evaluation metrics.
}

