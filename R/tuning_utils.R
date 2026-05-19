rank_hbic <- function(x, y, beta, n_original, df_thr = 0) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  beta <- as.numeric(beta)
  n_pairs <- length(y)
  p <- ncol(x)
  residual <- abs(as.numeric(y - x %*% beta))
  df <- sum(abs(beta) > df_thr)
  log(sum(residual) / n_pairs) +
    df * log(log(n_original)) * log(p) / n_original
}

.select_hbic_tol <- function(lambda_grid, hbic, tol_rel = 0.1,
                             prefer = c("small", "mid", "large")) {
  prefer <- match.arg(prefer)
  ord <- order(lambda_grid, decreasing = TRUE)
  lambda_grid <- lambda_grid[ord]
  hbic <- hbic[ord]
  id_min <- which.min(hbic)
  best <- hbic[id_min]
  threshold <- best + tol_rel * abs(best)
  candidates <- which(hbic <= threshold)
  id <- switch(
    prefer,
    large = min(candidates),
    small = max(candidates),
    mid = {
      mid <- (length(lambda_grid) + 1) / 2
      candidates[which.min(abs(candidates - mid))]
    }
  )
  list(id = id, lambda = lambda_grid[id], id_min = id_min,
       threshold = threshold)
}

rank_admm_select_lambda <- function(
    x, y,
    lambda_grid,
    penalty = c("L1", "SCAD"),
    M = 2L,
    G = 4L,
    phi = 1 / 400,
    rho = 1.618,
    max_iter = 500L,
    e.abs = 1e-2,
    e.rel = 1e-2,
    stop_rule = "residual",
    beta_init = NULL,
    rank_transform = TRUE,
    standardize = TRUE,
    n_original = NULL,
    scad_a = 3.7,
    lla_max = 2L,
    lla_tol = 1e-3,
    lla_damp = 0.25,
    lla_eps = 1e-6,
    l1_init_lambda_grid = NULL,
    scad_tol_rel = 0.1,
    scad_prefer = "small",
    ...
) {
  penalty <- match.arg(penalty)
  prep <- prepare_rank_data(x, y, standardize = standardize,
                            rank_transform = rank_transform)
  if (is.null(n_original)) n_original <- prep$n_original
  lambda_grid <- sort(unique(lambda_grid), decreasing = TRUE)
  
  if (penalty == "SCAD" && is.null(beta_init) &&
      !is.null(l1_init_lambda_grid)) {
    pilot <- rank_admm_select_lambda(
      prep$x, prep$y, lambda_grid = l1_init_lambda_grid,
      penalty = "L1", M = M, G = G, phi = phi, rho = rho,
      max_iter = max_iter, e.abs = e.abs, e.rel = e.rel,
      stop_rule = "residual", beta_init = NULL,
      rank_transform = FALSE, standardize = FALSE,
      n_original = n_original, ...
    )
    beta_init <- pilot$best_model$beta
  }
  
  fits <- vector("list", length(lambda_grid))
  hbic <- numeric(length(lambda_grid))
  beta_current <- beta_init
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    if (penalty == "L1") {
      fits[[i]] <- rank_admm_fit(
        prep$x, prep$y, penalty = "L1", lambda = lam,
        M = M, G = G, phi = phi, rho = rho, max_iter = max_iter,
        e.abs = e.abs, e.rel = e.rel, stop_rule = stop_rule,
        beta_start = beta_current, rank_transform = FALSE,
        standardize = FALSE, ...
      )
      beta_current <- fits[[i]]$beta
    } else {
      fits[[i]] <- rank_admm_fit(
        prep$x, prep$y, penalty = "SCAD", lambda = lam,
        M = M, G = G, phi = phi, rho = rho, max_iter = max_iter,
        e.abs = e.abs, e.rel = e.rel, stop_rule = stop_rule,
        beta_start = beta_init, rank_transform = FALSE,
        standardize = FALSE, scad_a = scad_a, lla_max = lla_max,
        lla_tol = lla_tol, lla_damp = lla_damp, lla_eps = lla_eps,
        ...
      )
    }
    hbic[i] <- rank_hbic(prep$x, prep$y, fits[[i]]$beta,
                         n_original = n_original)
  }
  
  if (penalty == "SCAD") {
    sel <- .select_hbic_tol(lambda_grid, hbic, tol_rel = scad_tol_rel,
                            prefer = scad_prefer)
    id <- sel$id
  } else {
    id <- which.min(hbic)
    sel <- NULL
  }
  list(best_lambda = lambda_grid[id], best_model = fits[[id]],
       hbic = hbic, lambda_grid = lambda_grid, id = id,
       selection = sel, n_original = n_original,
       x = prep$x, y = prep$y, x_original = prep$x_original,
       y_original = prep$y_original)
}
