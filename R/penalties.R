.soft_threshold <- function(x, lambda) {
  sign(x) * pmax(0, abs(x) - lambda)
}

.l2norm <- function(x) sqrt(sum(x^2))

.scad_penalty_value <- function(z, lambda, a = 3.7) {
  t <- abs(z)
  out <- numeric(length(t))
  idx1 <- t <= lambda
  idx2 <- (t > lambda) & (t <= a * lambda)
  idx3 <- t > a * lambda
  out[idx1] <- lambda * t[idx1]
  out[idx2] <- (-t[idx2]^2 + 2 * a * lambda * t[idx2] - lambda^2) /
    (2 * (a - 1))
  out[idx3] <- ((a + 1) * lambda^2) / 2
  out
}

.scad_deriv <- function(t, lambda, a = 3.7) {
  t <- pmax(t, 0)
  out <- numeric(length(t))
  out[t <= lambda] <- lambda
  idx <- (t > lambda) & (t <= a * lambda)
  out[idx] <- (a * lambda - t[idx]) / (a - 1)
  out[t > a * lambda] <- 0
  out
}

.scad_weights <- function(beta, lambda, a = 3.7, eps = 1e-6,
                          penalize_intercept = TRUE) {
  w <- .scad_deriv(abs(beta) + eps, lambda = lambda, a = a)
  if (!penalize_intercept && length(w) >= 1L) w[1L] <- 0
  w
}

.support_set <- function(beta, delta = 1e-6,
                         penalize_intercept = TRUE) {
  idx <- which(abs(beta) > delta)
  if (!penalize_intercept && length(idx) > 0L) idx <- setdiff(idx, 1L)
  idx
}
