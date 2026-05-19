.matsplit_mg <- function(x, M, G) {
  n <- nrow(x)
  p <- ncol(x)
  if (n %% M != 0L) stop("nrow(x) must be divisible by M", call. = FALSE)
  if (p %% G != 0L) stop("ncol(x) must be divisible by G", call. = FALSE)
  nr <- n %/% M
  nc <- p %/% G
  out <- matrix(vector("list", M * G), nrow = M, ncol = G)
  for (m in seq_len(M)) {
    r_idx <- ((m - 1L) * nr + 1L):(m * nr)
    for (g in seq_len(G)) {
      c_idx <- ((g - 1L) * nc + 1L):(g * nc)
      out[[m, g]] <- x[r_idx, c_idx, drop = FALSE]
    }
  }
  out
}

.rank_admm_l1_core <- function(
    x, y, lambda,
    beta_start = NULL,
    max_iter = 500L,
    rho = 1.618,
    tau = 0.5,
    e.rel = 1e-2,
    e.abs = 1e-2,
    phi = 1 / 400,
    M = 2L,
    G = 4L,
    K = 1L,
    stop_rule = c("residual", "objective", "beta", "support", "hybrid"),
    epsF = 1e-6,
    epsBeta = 1e-4,
    delta_support = 1e-6,
    L_support = 5L,
    penalize_intercept = TRUE,
    min_iter = 1L,
    residual_gated = TRUE,
    obj_scale = c("mean", "sum"),
    freeze_phi = TRUE,
    verbose = FALSE
) {
  stop_rule <- match.arg(stop_rule)
  obj_scale <- match.arg(obj_scale)
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)
  p <- ncol(x)
  if (length(y) != n) stop("length(y) must match nrow(x)", call. = FALSE)
  if (!(length(lambda) == 1L || length(lambda) == p)) {
    stop("lambda must be scalar or length ncol(x)", call. = FALSE)
  }
  if (any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("lambda must be finite and nonnegative", call. = FALSE)
  }
  if (is.null(beta_start) || length(beta_start) != p) beta_start <- numeric(p)
  
  x.mg <- .matsplit_mg(x, M = M, G = G)
  nr <- n %/% M
  nc <- p %/% G
  
  y.m <- matrix(0, nrow = nr, ncol = M)
  for (m in seq_len(M)) {
    r_idx <- ((m - 1L) * nr + 1L):(m * nr)
    y.m[, m] <- y[r_idx]
  }
  
  beta.g <- matrix(beta_start, nrow = nc, ncol = G)
  beta.mg <- vector("list", G)
  for (g in seq_len(G)) {
    beta.mg[[g]] <- matrix(0, nrow = nc, ncol = M)
    for (m in seq_len(M)) beta.mg[[g]][, m] <- beta.g[, g]
  }
  
  w.mg <- vector("list", M)
  z.m <- matrix(0, nrow = nr, ncol = M)
  for (m in seq_len(M)) {
    w.mg[[m]] <- matrix(0, nrow = nr, ncol = G)
    for (g in seq_len(G)) {
      w.mg[[m]][, g] <- as.numeric(x.mg[[m, g]] %*% beta.mg[[g]][, m])
    }
    z.m[, m] <- y.m[, m] - rowSums(w.mg[[m]])
  }
  
  theta.m <- matrix(0, nrow = nr, ncol = M)
  eta.mg <- vector("list", M)
  for (m in seq_len(M)) eta.mg[[m]] <- matrix(0, nrow = nr, ncol = G)
  gamma.g <- vector("list", G)
  for (g in seq_len(G)) gamma.g[[g]] <- matrix(0, nrow = nc, ncol = M)
  
  x.inv <- matrix(vector("list", M * G), nrow = M, ncol = G)
  for (m in seq_len(M)) for (g in seq_len(G)) {
    Xm <- x.mg[[m, g]]
    x.inv[[m, g]] <- solve(crossprod(Xm) + diag(ncol(Xm)))
  }
  
  beta_vec_prev <- c(beta.g)
  objective <- function(beta_vec) {
    r <- as.numeric(y - x %*% beta_vec)
    loss <- if (obj_scale == "mean") mean(abs(r)) else sum(abs(r))
    pen_idx <- seq_along(beta_vec)
    if (!penalize_intercept && length(pen_idx) >= 1L) pen_idx <- pen_idx[-1L]
    pen <- if (length(lambda) == 1L) {
      lambda * sum(abs(beta_vec[pen_idx]))
    } else {
      sum(lambda[pen_idx] * abs(beta_vec[pen_idx]))
    }
    loss + pen
  }
  F_prev <- objective(beta_vec_prev)
  A_prev <- .support_set(beta_vec_prev, delta = delta_support,
                         penalize_intercept = penalize_intercept)
  stable_ct <- 0L
  stop_iter <- NA_integer_
  stop_reason <- NA_character_
  
  prox_abs_phi <- function(v, phi) {
    thr <- if (obj_scale == "mean") 1 / (phi * nrow(x)) else 1 / phi
    .soft_threshold(v, thr)
  }
  lambda_slice_g <- function(g) {
    idx <- ((g - 1L) * nc + 1L):(g * nc)
    lam_g <- if (length(lambda) == 1L) rep(lambda, nc) else lambda[idx]
    if (!penalize_intercept) lam_g[idx == 1L] <- 0
    lam_g
  }
  
  t <- 1L
  r.norm <- s.norm <- C_pri <- C_dual <- delta_F <- delta_beta <- NA_real_
  while (t <= max_iter) {
    z_prev <- c(z.m)
    
    for (m in seq_len(M)) for (g in seq_len(G)) {
      Xm <- x.mg[[m, g]]
      rhs <- crossprod(Xm, w.mg[[m]][, g] - eta.mg[[m]][, g] / phi) +
        beta.g[, g] + gamma.g[[g]][, m] / phi
      beta.mg[[g]][, m] <- as.numeric(x.inv[[m, g]] %*% rhs)
    }
    
    for (g in seq_len(G)) {
      u_vec <- rowMeans(beta.mg[[g]]) - rowMeans(gamma.g[[g]]) / phi
      beta.g[, g] <- .soft_threshold(u_vec, lambda_slice_g(g) / (phi * M))
    }
    
    for (m in seq_len(M)) {
      v <- y.m[, m] - rowSums(w.mg[[m]]) - theta.m[, m] / phi
      z.m[, m] <- prox_abs_phi(v, phi)
    }
    
    for (m in seq_len(M)) {
      A_mg <- matrix(0, nrow = nr, ncol = G)
      for (g in seq_len(G)) {
        A_mg[, g] <- as.numeric(x.mg[[m, g]] %*% beta.mg[[g]][, m]) +
          eta.mg[[m]][, g] / phi
      }
      B_m <- y.m[, m] - z.m[, m] - theta.m[, m] / phi
      common_term <- (rowSums(A_mg) - B_m) / (1 + G)
      for (g in seq_len(G)) w.mg[[m]][, g] <- A_mg[, g] - common_term
    }
    
    for (m in seq_len(M)) {
      theta.m[, m] <- theta.m[, m] +
        rho * phi * (z.m[, m] + rowSums(w.mg[[m]]) - y.m[, m])
    }
    for (m in seq_len(M)) for (g in seq_len(G)) {
      Xm <- x.mg[[m, g]]
      eta.mg[[m]][, g] <- eta.mg[[m]][, g] +
        rho * phi * (as.numeric(Xm %*% beta.mg[[g]][, m]) - w.mg[[m]][, g])
    }
    for (g in seq_len(G)) for (m in seq_len(M)) {
      gamma.g[[g]][, m] <- gamma.g[[g]][, m] +
        rho * phi * (beta.g[, g] - beta.mg[[g]][, m])
    }
    
    r1 <- 0
    for (g in seq_len(G)) for (m in seq_len(M)) {
      r1 <- r1 + .l2norm(beta.mg[[g]][, m] - beta.g[, g])
    }
    r2 <- 0
    for (m in seq_len(M)) {
      r2 <- r2 + .l2norm(z.m[, m] + rowSums(w.mg[[m]]) - y.m[, m])
    }
    r3 <- 0
    for (m in seq_len(M)) for (g in seq_len(G)) {
      Xm <- x.mg[[m, g]]
      r3 <- r3 + .l2norm(as.numeric(Xm %*% beta.mg[[g]][, m]) -
                           w.mg[[m]][, g])
    }
    r.norm <- r1 + r2 + r3
    s.norm <- .l2norm(phi * crossprod(x, c(z.m) - z_prev))
    C_pri <- sqrt(n) * e.abs +
      e.rel * (max(.l2norm(x %*% c(beta.g)), .l2norm(c(z.m)),
                   .l2norm(y)) / sqrt(M * G))
    C_dual <- sqrt(p) * e.abs +
      e.rel * (.l2norm(crossprod(x, c(theta.m))) / (M * G))
    
    if (!freeze_phi && (t %% K) == 0L && t < 100L) {
      phi_old <- phi
      if (r.norm > 3 * s.norm) phi <- 2 * phi
      else if (s.norm > 3 * r.norm) phi <- phi / 2
      if (!identical(phi, phi_old)) {
        sc <- phi_old / phi
        for (g in seq_len(G)) gamma.g[[g]] <- gamma.g[[g]] * sc
        theta.m <- theta.m * sc
        for (m in seq_len(M)) eta.mg[[m]] <- eta.mg[[m]] * sc
      }
    }
    
    beta_vec <- c(beta.g)
    F_curr <- objective(beta_vec)
    delta_beta <- max(abs(beta_vec - beta_vec_prev)) /
      max(1, max(abs(beta_vec_prev)))
    delta_F <- abs(F_curr - F_prev) / max(1, abs(F_prev))
    A_curr <- .support_set(beta_vec, delta = delta_support,
                           penalize_intercept = penalize_intercept)
    stable_ct <- if (identical(A_curr, A_prev)) stable_ct + 1L else 0L
    feas_ok <- isTRUE(r.norm <= C_pri && s.norm <= C_dual)
    gate <- t >= min_iter
    stop_now <- switch(
      stop_rule,
      residual = gate && feas_ok,
      objective = {
        core <- delta_F <= epsF
        if (residual_gated) gate && feas_ok && core else gate && core
      },
      beta = {
        core <- delta_beta <= epsBeta
        if (residual_gated) gate && feas_ok && core else gate && core
      },
      support = {
        core <- stable_ct >= L_support
        if (residual_gated) gate && feas_ok && core else gate && core
      },
      hybrid = gate && feas_ok && delta_beta <= epsBeta
    )
    if (isTRUE(stop_now)) {
      stop_iter <- t
      stop_reason <- stop_rule
      break
    }
    beta_vec_prev <- beta_vec
    F_prev <- F_curr
    A_prev <- A_curr
    t <- t + 1L
  }
  if (is.na(stop_iter)) {
    stop_iter <- t
    stop_reason <- "max_iter"
  }
  structure(
    list(beta = c(beta.g), beta.hat = c(beta.g), lambda = lambda,
         iter = t, stop_iter = stop_iter, stop_rule = stop_rule,
         stop_reason = stop_reason, converged = isTRUE(r.norm <= C_pri &&
                                                         s.norm <= C_dual),
         feas_ok = isTRUE(r.norm <= C_pri && s.norm <= C_dual),
         r_norm = r.norm, s_norm = s.norm, e_pri = C_pri, e_dual = C_dual,
         delta_F = delta_F, delta_beta = delta_beta, phi_final = phi),
    class = "rank_admm_fit"
  )
}

.fit_scad_lla <- function(x, y, lambda, beta_init = NULL, a = 3.7,
                          lla_max = 2L, lla_tol = 1e-3,
                          lla_damp = 0.25, eps = 1e-6,
                          penalize_intercept = TRUE,
                          monotone = TRUE,
                          line_search_min = 1 / 64,
                          line_search_shrink = 0.5,
                          obj_tol = 1e-10,
                          ...) {
  p <- ncol(x)
  if (is.null(beta_init) || length(beta_init) != p) beta_init <- numeric(p)
  if (!is.numeric(lla_damp) || length(lla_damp) != 1L ||
      !is.finite(lla_damp) || lla_damp <= 0 || lla_damp > 1) {
    stop("lla_damp must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(line_search_min) || length(line_search_min) != 1L ||
      !is.finite(line_search_min) || line_search_min <= 0 ||
      line_search_min > 1) {
    stop("line_search_min must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(line_search_shrink) ||
      length(line_search_shrink) != 1L ||
      !is.finite(line_search_shrink) || line_search_shrink <= 0 ||
      line_search_shrink >= 1) {
    stop("line_search_shrink must be in (0, 1).", call. = FALSE)
  }

  dots <- list(...)
  obj_scale <- if (!is.null(dots$obj_scale)) as.character(dots$obj_scale) else "mean"
  obj_scale <- match.arg(obj_scale, c("mean", "sum"))
  scad_objective <- function(beta_vec) {
    r <- as.numeric(y - x %*% beta_vec)
    loss <- if (obj_scale == "mean") mean(abs(r)) else sum(abs(r))
    pen_idx <- seq_along(beta_vec)
    if (!penalize_intercept && length(pen_idx) >= 1L) pen_idx <- pen_idx[-1L]
    loss + sum(.scad_penalty_value(beta_vec[pen_idx], lambda = lambda, a = a))
  }

  beta_old <- beta_init
  obj_old <- scad_objective(beta_old)
  fit <- NULL
  for (it in seq_len(lla_max)) {
    w <- .scad_weights(beta_old, lambda = lambda, a = a, eps = eps,
                       penalize_intercept = penalize_intercept)
    fit <- .rank_admm_l1_core(
      x = x, y = y, lambda = w, beta_start = beta_old,
      penalize_intercept = penalize_intercept, ...
    )
    beta_candidate <- as.numeric(fit$beta)
    step <- lla_damp
    beta_new <- (1 - step) * beta_old + step * beta_candidate
    obj_new <- scad_objective(beta_new)

    ## Equation-level link: LLA majorizes the SCAD penalty by a weighted L1
    ## problem. With loose inner ADMM solves, the line search keeps the actual
    ## SCAD rank objective monotone without changing the estimator definition.
    if (isTRUE(monotone)) {
      while ((!is.finite(obj_new) || obj_new > obj_old + obj_tol) &&
             step > line_search_min) {
        step <- step * line_search_shrink
        beta_new <- (1 - step) * beta_old + step * beta_candidate
        obj_new <- scad_objective(beta_new)
      }
      if (!is.finite(obj_new) || obj_new > obj_old + obj_tol) {
        beta_new <- beta_old
        obj_new <- obj_old
        step <- 0
      }
    }

    relchg <- max(abs(beta_new - beta_old)) / max(1, max(abs(beta_old)))
    beta_old <- beta_new
    obj_old <- obj_new
    fit$lla_step <- step
    fit$scad_objective <- obj_new
    if (is.finite(relchg) && relchg < lla_tol) break
  }
  fit$beta <- beta_old
  fit$beta.hat <- beta_old
  fit$lambda <- lambda
  fit$lla_iter <- it
  fit$lla_damp <- lla_damp
  fit$penalty <- "SCAD"
  class(fit) <- "rank_admm_fit"
  fit
}

rank_admm_fit <- function(
    x, y,
    penalty = c("L1", "SCAD"),
    lambda,
    M = 2L,
    G = 4L,
    phi = 1 / 400,
    rho = 1.618,
    max_iter = 500L,
    e.abs = 1e-2,
    e.rel = 1e-2,
    stop_rule = c("residual", "objective", "beta", "support", "hybrid"),
    beta_start = NULL,
    rank_transform = FALSE,
    standardize = FALSE,
    tau = 0.5,
    K = 1L,
    obj_scale = c("mean", "sum"),
    freeze_phi = TRUE,
    penalize_intercept = TRUE,
    scad_a = 3.7,
    lla_max = 2L,
    lla_tol = 1e-3,
    lla_damp = 0.25,
    lla_eps = 1e-6,
    scad_monotone = TRUE,
    scad_line_search_min = 1 / 64,
    scad_line_search_shrink = 0.5,
    scad_obj_tol = 1e-10,
    ...
) {
  penalty <- match.arg(penalty)
  stop_rule <- match.arg(stop_rule)
  obj_scale <- match.arg(obj_scale)
  prep <- prepare_rank_data(x, y, standardize = standardize,
                            rank_transform = rank_transform)
  fit_args <- list(
    x = prep$x, y = prep$y, lambda = lambda, beta_start = beta_start,
    max_iter = max_iter, rho = rho, tau = tau, e.rel = e.rel,
    e.abs = e.abs, phi = phi, M = M, G = G, K = K,
    stop_rule = stop_rule, obj_scale = obj_scale, freeze_phi = freeze_phi,
    penalize_intercept = penalize_intercept, ...
  )
  fit <- if (penalty == "L1") {
    do.call(.rank_admm_l1_core, fit_args)
  } else {
    fit_args$beta_init <- fit_args$beta_start
    fit_args$beta_start <- NULL
    fit_args$a <- scad_a
    fit_args$lla_max <- lla_max
    fit_args$lla_tol <- lla_tol
    fit_args$lla_damp <- lla_damp
    fit_args$eps <- lla_eps
    fit_args$monotone <- scad_monotone
    fit_args$line_search_min <- scad_line_search_min
    fit_args$line_search_shrink <- scad_line_search_shrink
    fit_args$obj_tol <- scad_obj_tol
    do.call(.fit_scad_lla, fit_args)
  }
  fit$penalty <- penalty
  fit$n_original <- prep$n_original
  fit$x_original <- prep$x_original
  fit$y_original <- prep$y_original
  fit
}
