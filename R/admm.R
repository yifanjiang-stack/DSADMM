# DS-ADMM solver with multiple stopping rules

admm_both_stop <- function(
    X, y,
    lambda = 0,
    max_iter = 5000,
    rho = 1.618,
    tau = 0.5,
    e.rel = 1e-3,
    e.abs = 1e-3,
    phi = 1.0,
    M = 2,
    G = 1,
    K = 1,
    stop_rule = c("residual", "objective", "beta", "support", "hybrid"),
    epsF = 1e-6,
    epsBeta = 1e-4,
    delta_support = 1e-6,
    L_support = 5,
    penalize_intercept = FALSE,
    min_iter = 10,
    residual_gated = TRUE,
    obj_scale = c("mean", "sum"),
    freeze_phi = FALSE,
    verbose = FALSE,
    print_every = 50,
    beta_start = NULL,
    return_trace = FALSE,
    trace_every = 1,
    trace_max = 2000
) {
  stop_rule <- match.arg(stop_rule)
  obj_scale <- match.arg(obj_scale)

  stopifnot(is.matrix(X), is.numeric(y), length(y) == nrow(X))
  n <- nrow(X); p <- ncol(X)
  stopifnot(n %% M == 0, p %% G == 0)

  # Scale threshold by number of rows (N)
  N_pairs_scale <- n

  x.mg <- matsplit_mg(X, M, G)

  nr <- n %/% M
  nc <- p %/% G

  ## y block split
  y.m <- matrix(0, nrow = nr, ncol = M)
  for (m in 1:M) {
    r_idx <- ((m - 1) * nr + 1):(m * nr)
    y.m[, m] <- y[r_idx]
  }

  # Initialization
  if (is.null(beta_start) || length(beta_start) != p) beta_start <- numeric(p)

  beta.g <- matrix(beta_start, nrow = nc, ncol = G)

  beta.mg <- vector("list", G)
  for (g in 1:G) {
    beta.mg[[g]] <- matrix(0, nrow = nc, ncol = M)
    for (m in 1:M) beta.mg[[g]][, m] <- beta.g[, g]
  }

  w.mg <- vector("list", M)
  z.m <- matrix(0, nrow = nr, ncol = M)

  for (m in 1:M) {
    w.mg[[m]] <- matrix(0, nrow = nr, ncol = G)
    for (g in 1:G) {
      w.mg[[m]][, g] <- as.numeric(x.mg[[m, g]] %*% beta.mg[[g]][, m])
    }
    z.m[, m] <- y.m[, m] - rowSums(w.mg[[m]])
  }

  ## Dual variables
  theta.m <- matrix(0, nrow = nr, ncol = M)
  eta.mg <- vector("list", M)
  for (m in 1:M) eta.mg[[m]] <- matrix(0, nrow = nr, ncol = G)
  gamma.g <- vector("list", G)
  for (g in 1:G) gamma.g[[g]] <- matrix(0, nrow = nc, ncol = M)

  ## Precompute inverse
  x.inv <- matrix(vector("list", M * G), nrow = M, ncol = G)
  for (m in 1:M) for (g in 1:G) {
    Xm <- x.mg[[m, g]]
    x.inv[[m, g]] <- solve(crossprod(Xm) + diag(ncol(Xm)))
  }

  beta_vec_prev <- c(beta.g)

  obj_lad_l1 <- function(beta_vec) {
    r <- as.numeric(y - X %*% beta_vec)
    loss <- sum(abs(r))
    pen_idx <- seq_along(beta_vec)
    if (!penalize_intercept && length(beta_vec) >= 1) pen_idx <- pen_idx[-1]
    pen <- lambda * sum(abs(beta_vec[pen_idx]))
    loss + pen
  }

  F_prev <- obj_lad_l1(beta_vec_prev)
  if (!is.finite(F_prev)) F_prev <- 0

  A_prev <- support_set(beta_vec_prev, delta = delta_support, penalize_intercept = penalize_intercept)
  stable_ct <- 0L
  stop_iter <- NA_integer_
  stop_reason <- NA_character_

  # trace placeholder
  trace <- NULL

  prox_abs_phi <- function(v, phi, scale_factor) {
    thr <- 1 / (phi * scale_factor)
    sign(v) * pmax(abs(v) - thr, 0)
  }

  ## Initialize loop-updated quantities so the return list is always defined.
  r.norm <- NA_real_
  s.norm <- NA_real_
  C_pri  <- NA_real_
  C_dual <- NA_real_
  delta_F <- NA_real_
  delta_beta <- NA_real_

  t <- 1L
  while (t < max_iter) {
    z_prev <- c(z.m)

    ## (1) beta.mg update
    for (m in 1:M) for (g in 1:G) {
      Xm <- x.mg[[m, g]]
      rhs <- crossprod(Xm, w.mg[[m]][, g] - eta.mg[[m]][, g] / phi) +
        beta.g[, g] + gamma.g[[g]][, m] / phi
      beta.mg[[g]][, m] <- as.numeric(x.inv[[m, g]] %*% rhs)
    }

    ## (2) beta.g update
    for (g in 1:G) {
      u_vec <- rowMeans(beta.mg[[g]]) - (1 / phi) * rowMeans(gamma.g[[g]])
      a_val <- lambda / (phi * M)
      beta.g[, g] <- prox_beta(u_vec, a_val)
    }

    ## (3) z.m update
    for (m in 1:M) {
      v <- y.m[, m] - rowSums(w.mg[[m]]) - theta.m[, m] / phi
      z.m[, m] <- prox_abs_phi(v, phi, scale_factor = N_pairs_scale)
    }

    ## (4) w.mg update (exact with duals)
    for (m in 1:M) {
      A_mg <- matrix(0, nrow = nr, ncol = G)
      for (g in 1:G) {
        A_mg[, g] <- as.numeric(x.mg[[m, g]] %*% beta.mg[[g]][, m]) +
          eta.mg[[m]][, g] / phi
      }

      B_m <- y.m[, m] - z.m[, m] - theta.m[, m] / phi
      sum_A <- rowSums(A_mg)
      common_term <- (sum_A - B_m) / (1 + G)

      for (g in 1:G) {
        w.mg[[m]][, g] <- A_mg[, g] - common_term
      }
    }

    ## (5) Multipliers update
    for (m in 1:M) {
      theta.m[, m] <- theta.m[, m] + rho * phi * (z.m[, m] + rowSums(w.mg[[m]]) - y.m[, m])
    }
    for (m in 1:M) for (g in 1:G) {
      Xm <- x.mg[[m, g]]
      eta.mg[[m]][, g] <- eta.mg[[m]][, g] + rho * phi * (as.numeric(Xm %*% beta.mg[[g]][, m]) - w.mg[[m]][, g])
    }
    for (g in 1:G) for (m in 1:M) {
      gamma.g[[g]][, m] <- gamma.g[[g]][, m] + rho * phi * (beta.g[, g] - beta.mg[[g]][, m])
    }

    ## (6) Residuals
    r1 <- 0
    for (g in 1:G) for (m in 1:M) r1 <- r1 + l2norm(beta.mg[[g]][, m] - beta.g[, g])
    r2 <- 0
    for (m in 1:M) r2 <- r2 + l2norm(z.m[, m] + rowSums(w.mg[[m]]) - y.m[, m])
    r3 <- 0
    for (m in 1:M) for (g in 1:G) {
      Xm <- x.mg[[m, g]]
      r3 <- r3 + l2norm(as.numeric(Xm %*% beta.mg[[g]][, m]) - w.mg[[m]][, g])
    }
    r.norm <- r1 + r2 + r3
    s <- phi * crossprod(X, (c(z.m) - z_prev))
    s.norm <- l2norm(s)

    ## Thresholds
    C_pri <- sqrt(n) * e.abs + e.rel * (max(l2norm(X %*% c(beta.g)), l2norm(c(z.m)), l2norm(y)) / sqrt(M * G))
    C_dual <- sqrt(p) * e.abs + e.rel * (l2norm(crossprod(X, c(theta.m))) / (M * G))

    ## Adaptive phi
    if (!freeze_phi && (t %% K) == 0 && t < 100) {
      phi_old <- phi
      if (r.norm > 10 * s.norm) {
        phi <- 10 * phi
      } else if (s.norm > 10 * r.norm) {
        phi <- phi / 10
      }
      if (!identical(phi, phi_old)) {
        sc <- phi_old / phi
        for (g in 1:G) gamma.g[[g]] <- gamma.g[[g]] * sc
        theta.m <- theta.m * sc
        for (m in 1:M) eta.mg[[m]] <- eta.mg[[m]] * sc
      }
    }

    ## (7) Stopping checks
    beta_vec <- c(beta.g)
    F_curr <- obj_lad_l1(beta_vec)
    if (!is.finite(F_curr)) F_curr <- F_prev

    denom_b <- max(1, max(abs(beta_vec_prev)))
    delta_beta <- max(abs(beta_vec - beta_vec_prev)) / denom_b
    if (!is.finite(delta_beta)) delta_beta <- Inf

    delta_F <- abs(F_curr - F_prev) / max(1, abs(F_prev))
    if (!is.finite(delta_F)) delta_F <- Inf

    A_curr <- support_set(beta_vec, delta = delta_support, penalize_intercept = penalize_intercept)
    if (identical(A_curr, A_prev)) stable_ct <- stable_ct + 1L else stable_ct <- 0L

    feas_ok <- isTRUE(r.norm <= C_pri && s.norm <= C_dual)
    gate <- (t >= min_iter)

    stop_now <- switch(stop_rule,
      residual = gate && feas_ok,
      objective = {
        core <- (delta_F <= epsF)
        if (residual_gated) (gate && feas_ok && core) else (gate && core)
      },
      beta = {
        core <- (delta_beta <= epsBeta)
        if (residual_gated) (gate && feas_ok && core) else (gate && core)
      },
      support = {
        core <- (stable_ct >= L_support)
        if (residual_gated) (gate && feas_ok && core) else (gate && core)
      },
      hybrid = {
        core <- (delta_beta <= epsBeta)
        (gate && feas_ok && core)
      }
    )
    if (!isTRUE(stop_now)) stop_now <- FALSE

    if (stop_now) {
      stop_iter <- t
      stop_reason <- stop_rule
      break
    }

    beta_vec_prev <- beta_vec
    F_prev <- F_curr
    A_prev <- A_curr
    t <- t + 1L

    if (verbose && (t %% print_every == 0)) {
      message("iter=", t, " r_norm=", signif(r.norm, 4), " s_norm=", signif(s.norm, 4))
    }
  }

  if (is.na(stop_iter)) { stop_iter <- t; stop_reason <- "max_iter" }

  list(
    beta.hat = c(beta.g),
    iter = t,
    stop_iter = stop_iter,
    stop_rule = stop_rule,
    stop_reason = stop_reason,
    r_norm = r.norm,
    s_norm = s.norm,
    e_pri = C_pri,
    e_dual = C_dual,
    feas_ok = isTRUE(r.norm <= C_pri && s.norm <= C_dual),
    F = F_prev,
    delta_F = delta_F,
    delta_beta = delta_beta,
    support_size = length(A_prev),
    phi_final = phi,
    trace = trace
  )
}
