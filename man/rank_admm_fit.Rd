\name{rank_admm_fit}
\alias{rank_admm_fit}
\title{Fit Penalized Rank Regression by Double-Splitting ADMM}
\usage{
rank_admm_fit(x, y, penalty = c("L1", "SCAD"), lambda, M = 2L,
  G = 4L, phi = 1/400, rho = 1.618, max_iter = 500L,
  e.abs = 1e-2, e.rel = 1e-2, stop_rule = c("residual",
  "objective", "beta", "support", "hybrid"), beta_start = NULL,
  rank_transform = FALSE, standardize = FALSE, tau = 0.5, K = 1L,
  obj_scale = c("mean", "sum"), freeze_phi = TRUE,
  penalize_intercept = TRUE, scad_a = 3.7, lla_max = 2L,
  lla_tol = 1e-3, lla_damp = 0.25, lla_eps = 1e-6,
  scad_monotone = TRUE, scad_line_search_min = 1/64,
  scad_line_search_shrink = 0.5, scad_obj_tol = 1e-10, ...)
}
\arguments{
\item{x}{Design matrix. Use raw data with \code{rank_transform = TRUE}, or a pairwise-difference design with \code{rank_transform = FALSE}.}
\item{y}{Response vector matching \code{x}.}
\item{penalty}{Penalty type, either \code{"L1"} or \code{"SCAD"}.}
\item{lambda}{Penalty tuning parameter. For L1 this may be scalar or coefficient-specific.}
\item{M, G, phi, rho}{ADMM splitting and penalty settings used in the simulation scripts.}
\item{max_iter}{Maximum ADMM iterations.}
\item{e.abs, e.rel}{Absolute and relative ADMM stopping tolerances.}
\item{stop_rule}{Stopping rule.}
\item{beta_start}{Optional starting coefficient vector.}
\item{rank_transform}{Whether to apply the pairwise rank-regression transformation.}
\item{standardize}{Whether to center/scale \code{x} and center \code{y} before fitting.}
\item{tau, K, obj_scale, freeze_phi, penalize_intercept}{Controls preserved from the original ADMM scripts.}
\item{scad_a, lla_max, lla_tol, lla_damp, lla_eps}{SCAD local-linear-approximation controls.}
\item{scad_monotone, scad_line_search_min, scad_line_search_shrink, scad_obj_tol}{SCAD objective monotonicity line-search controls.}
\item{...}{Additional ADMM or SCAD LLA controls.}
}
\value{A list with coefficient vector \code{beta}, tuning and convergence information, and class \code{rank_admm_fit}.}
\description{
Fits the L1 or SCAD penalized rank-regression objective using the double-splitting ADMM updates from the rank-ADMM scripts. SCAD uses a local-linear-approximation/MM path with ADMM inner solves.
}
