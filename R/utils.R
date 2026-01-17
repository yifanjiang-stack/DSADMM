# Utilities used by DS-ADMM and examples

creating_x <- function(x) {
  n <- nrow(x); p <- ncol(x)
  N_pairs <- n * (n - 1) / 2
  x_tilde <- matrix(0, nrow = N_pairs, ncol = p)
  k <- 1L
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      x_tilde[k, ] <- x[i, ] - x[j, ]
      k <- k + 1L
    }
  }
  x_tilde
}

creating_y <- function(y) {
  n <- length(y)
  N_pairs <- n * (n - 1) / 2
  y_tilde <- numeric(N_pairs)
  k <- 1L
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      y_tilde[k] <- y[i] - y[j]
      k <- k + 1L
    }
  }
  y_tilde
}

matsplit_mg <- function(X, M, G) {
  stopifnot(is.matrix(X))
  n <- nrow(X); p <- ncol(X)
  stopifnot(n %% M == 0, p %% G == 0)
  nr <- n %/% M
  nc <- p %/% G

  out <- matrix(vector("list", M * G), nrow = M, ncol = G)
  for (m in 1:M) {
    r_idx <- ((m - 1) * nr + 1):(m * nr)
    for (g in 1:G) {
      c_idx <- ((g - 1) * nc + 1):(g * nc)
      out[[m, g]] <- X[r_idx, c_idx, drop = FALSE]
    }
  }
  out
}

prox_beta <- function(u, a) sign(u) * pmax(0, abs(u) - a)

l2norm <- function(x) sqrt(sum(x^2))

AR <- function(n, p, rho) {
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  X <- matrix(0, nrow = n, ncol = p)
  X[, 1] <- z[, 1]
  for (j in 2:p) X[, j] <- rho * X[, j - 1] + sqrt(1 - rho^2) * z[, j]
  X
}

support_set <- function(beta, delta = 1e-6, penalize_intercept = FALSE) {
  j <- which(abs(beta) > delta)
  if (!penalize_intercept && length(j) > 0) j <- setdiff(j, 1L)
  j
}

threshold_beta <- function(beta, thr) {
  beta * (abs(beta) > thr)
}

calc_fp_fn <- function(beta_hat, beta_true, thr = 0) {
  p <- length(beta_true)
  sel <- which(abs(beta_hat) > thr)
  idx_true <- which(beta_true != 0)
  fp <- sum(!(sel %in% idx_true)) / max(1, length(setdiff(1:p, idx_true)))
  fn <- sum(!(idx_true %in% sel)) / max(1, length(idx_true))
  list(sel = sel, fp = fp, fn = fn)
}
