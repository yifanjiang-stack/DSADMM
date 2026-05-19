\name{rank_admm_select_lambda}
\alias{rank_admm_select_lambda}
\title{Select Lambda by HBIC}
\usage{
rank_admm_select_lambda(x, y, lambda_grid, penalty = c("L1", "SCAD"),
  M = 2L, G = 4L, phi = 1/400, rho = 1.618,
  max_iter = 500L, e.abs = 1e-2, e.rel = 1e-2,
  stop_rule = "residual", beta_init = NULL,
  rank_transform = TRUE, standardize = TRUE, n_original = NULL,
  scad_a = 3.7, lla_max = 2L, lla_tol = 1e-3,
  lla_damp = 0.25, lla_eps = 1e-6,
  l1_init_lambda_grid = NULL, scad_tol_rel = 0.1,
  scad_prefer = "small", ...)
}
\arguments{
\item{x}{Raw design matrix or pairwise-difference design.}
\item{y}{Response vector.}
\item{lambda_grid}{Candidate lambda values.}
\item{penalty}{Penalty type.}
\item{M, G, phi, rho, max_iter, e.abs, e.rel, stop_rule}{ADMM settings.}
\item{beta_init}{Optional starting coefficient vector.}
\item{rank_transform, standardize}{Preprocessing flags.}
\item{n_original}{Original sample size used in the HBIC complexity term.}
\item{scad_a, lla_max, lla_tol, lla_damp, lla_eps}{SCAD local-linear-approximation controls.}
\item{l1_init_lambda_grid}{Optional L1 pilot grid for SCAD initialization.}
\item{scad_tol_rel, scad_prefer}{Tolerance-selection controls for the SCAD HBIC path.}
\item{...}{Additional controls passed to \code{rank_admm_fit}.}
}
\value{A list containing \code{best_lambda}, \code{best_model}, HBIC values, and transformed data used for tuning.}
\description{
Runs the HBIC lambda-selection workflow used by the Example 1 and Example 2 simulation scripts.
}
