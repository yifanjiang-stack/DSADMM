# HBIC utilities

HBIC <- function(lamb, X, y, x.ori, beta_start = NULL,
                 stop_rule = c("residual","objective","beta","support","hybrid"),
                 freeze_phi = FALSE, tau = 0.5, M = 2, G = 1, K = 1, phi = 1.0, rho = 1.618,
                 max_iter = 100, e.rel = 1e-2, e.abs = 1e-2, epsF = 1e-3, epsBeta = 1e-3,
                 delta_support = 1e-4, L_support = 5, penalize_intercept = FALSE,
                 residual_gated = TRUE, obj_scale = "mean") {

  stop_rule <- match.arg(stop_rule)
  n_pen <- nrow(x.ori)
  p <- ncol(X)

  t0 <- proc.time()
  fit <- admm_both_stop(
    X = X, y = y, lambda = lamb, beta_start = beta_start, max_iter = max_iter,
    rho = rho, tau = tau, e.rel = e.rel, e.abs = e.abs, phi = phi, M = M, G = G, K = K,
    stop_rule = stop_rule, epsF = epsF, epsBeta = epsBeta, delta_support = delta_support,
    L_support = L_support, penalize_intercept = penalize_intercept, min_iter = 1,
    residual_gated = residual_gated, obj_scale = obj_scale, freeze_phi = freeze_phi,
    verbose = FALSE
  )
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])

  residual <- abs(as.numeric(y - X %*% fit$beta.hat))
  df <- length(which(fit$beta.hat != 0))
  N_pairs <- length(y)

  hbic_val <- log(sum(residual) / N_pairs) + df * log(log(n_pen)) * log(p) / n_pen

  list(hbic = hbic_val, model = fit, time = elapsed)
}

choose_lambda_hbic <- function(lamb_list, X, y, x.ori, beta_init_guess = NULL, ...) {
  lamb_list <- sort(lamb_list, decreasing = TRUE)
  hb <- numeric(length(lamb_list))
  models <- vector("list", length(lamb_list))
  times <- numeric(length(lamb_list))
  current_beta <- beta_init_guess

  for (i in seq_along(lamb_list)) {
    res <- HBIC(lamb = lamb_list[i], X = X, y = y, x.ori = x.ori, beta_start = current_beta, ...)
    hb[i] <- res$hbic
    models[[i]] <- res$model
    times[i] <- res$time
    current_beta <- res$model$beta.hat
  }

  id <- which.min(hb)
  list(best_lambda = lamb_list[id], hbic = hb, id = id,
       best_model = models[[id]], best_time = times[id])
}
