make_rank_pairs <- function(x, y) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)
  if (length(y) != n) stop("length(y) must match nrow(x)", call. = FALSE)
  if (n < 2L) stop("rank regression requires at least two observations",
                   call. = FALSE)
  
  n_pairs <- n * (n - 1L) / 2L
  x_pair <- matrix(0, nrow = n_pairs, ncol = p)
  y_pair <- numeric(n_pairs)
  k <- 1L
  for (i in seq_len(n - 1L)) {
    rows <- (i + 1L):n
    len <- length(rows)
    idx <- k:(k + len - 1L)
    x_pair[idx, ] <- matrix(x[i, ], nrow = len, ncol = p,
                            byrow = TRUE) - x[rows, , drop = FALSE]
    y_pair[idx] <- y[i] - y[rows]
    k <- k + len
  }
  list(x = x_pair, y = y_pair, n_original = n)
}

prepare_rank_data <- function(x, y, standardize = TRUE,
                              rank_transform = TRUE) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  if (length(y) != nrow(x)) {
    stop("length(y) must match nrow(x)", call. = FALSE)
  }
  if (standardize) {
    x <- scale(x, center = TRUE, scale = TRUE)
    y <- as.numeric(y - mean(y))
  }
  if (rank_transform) {
    pair <- make_rank_pairs(x, y)
    return(list(x = pair$x, y = pair$y, x_original = x,
                y_original = y, n_original = pair$n_original))
  }
  list(x = x, y = y, x_original = x, y_original = y, n_original = nrow(x))
}

simulate_rank_example <- function(n = 60L, p = 200L, rho = 0.5,
                                  error = c("normal", "t4", "mixture"),
                                  beta = NULL) {
  error <- match.arg(error)
  if (is.null(beta)) {
    beta <- numeric(p)
    beta[seq_len(min(4L, p))] <- 1
  }
  idx <- 0:(p - 1L)
  sigma <- rho^abs(outer(idx, idx, "-"))
  x <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma)
  eps <- switch(
    error,
    normal = stats::rnorm(n),
    t4 = sqrt(2) * stats::rt(n, df = 4),
    mixture = stats::rnorm(n, sd = stats::runif(n, min = 1, max = 5))
  )
  y <- as.numeric(x %*% beta + eps)
  list(x = x, y = y, beta = beta, error = error, rho = rho)
}

rank_admm_metrics <- function(beta_hat, beta_true, threshold = 0) {
  beta_hat <- as.numeric(beta_hat)
  beta_true <- as.numeric(beta_true)
  if (length(beta_hat) != length(beta_true)) {
    stop("beta_hat and beta_true must have the same length", call. = FALSE)
  }
  selected <- which(abs(beta_hat) > threshold)
  active <- which(beta_true != 0)
  inactive <- which(beta_true == 0)
  data.frame(
    size = length(selected),
    err_l1 = sum(abs(beta_hat - beta_true)),
    FP = sum(selected %in% inactive) / max(1L, length(inactive)),
    FN = sum(!(active %in% selected)) / max(1L, length(active))
  )
}
