## Example 1 / Table 1: L1-penalized rank ADMM.
## This script mirrors the updated simulation workflow:
## select lambda on the 500-iteration HBIC path, then report the same
## penalized L1 ADMM estimator at final max_iter = 400.

suppressPackageStartupMessages(library(DSADMM))

## PDF/script simulation setting.
n0 <- 60L
p0 <- 200L
rho0 <- 0.5
sigma0 <- 1
beta0 <- numeric(p0)
beta0[1:4] <- 1

## Original ADMM/tuning settings preserved from Example 1.R.
lamb_list_resid <- seq(0.15, 0.5, by = 0.01)
M0 <- 2L
G0 <- 4L
K0 <- 1L
phi0 <- 1 / 400
rho_admm <- 1.618
HBIC_MAXITER_res <- 500L
FINAL_MAXITER_res <- 400L
HBIC_EABS <- 1e-2
OBJ_SCALE <- "mean"

n_rep <- as.integer(Sys.getenv("DSADMM_EXAMPLE_REPS", "1"))
seed_master <- 2021L
seeds <- seed_master + seq_len(n_rep)
run_competitors <- !identical(Sys.getenv("DSADMM_RUN_COMPETITORS"), "FALSE")

glmnet_start <- function(x, y) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    if (identical(Sys.getenv("DSADMM_ALLOW_ZERO_START"), "TRUE")) {
      warning("glmnet is unavailable; using zero start for this smoke run.")
      return(numeric(ncol(x)))
    }
    stop("Example 1 needs glmnet to preserve the original cold-start rule. ",
         "Install glmnet or set DSADMM_ALLOW_ZERO_START=TRUE for a reduced smoke run.",
         call. = FALSE)
  }
  cv <- glmnet::cv.glmnet(x, y, alpha = 1, intercept = TRUE)
  as.numeric(stats::coef(cv, s = "lambda.min"))[-1]
}

make_row <- function(rep_id, seed, method, lambda, beta_hat) {
  met <- rank_admm_metrics(beta_hat, beta0, threshold = 0)
  data.frame(
    rep = rep_id,
    seed = seed,
    method = method,
    lambda = lambda,
    size = met$size,
    err_l1 = met$err_l1,
    FP = met$FP,
    FN = met$FN,
    V1 = beta_hat[1], V2 = beta_hat[2],
    V3 = beta_hat[3], V4 = beta_hat[4],
    stringsAsFactors = FALSE
  )
}

na_row <- function(rep_id, seed, method) {
  data.frame(
    rep = rep_id, seed = seed, method = method, lambda = NA_real_,
    size = NA_integer_, err_l1 = NA_real_, FP = NA_real_, FN = NA_real_,
    V1 = NA_real_, V2 = NA_real_, V3 = NA_real_, V4 = NA_real_,
    stringsAsFactors = FALSE
  )
}

safe_run <- function(expr, rep_id, seed, method) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[SKIP] %s failed at rep=%s seed=%s: %s",
                    method, rep_id, seed, e$message))
    na_row(rep_id, seed, method)
  })
}

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

check_loss <- function(r, tau) sum(r * (tau - (r < 0)))

shrinkcpp1 <- function(u, v) {
  (1 + sign(u - v)) / 2 * (u - v) -
    (1 + sign(-u - v)) / 2 * (-u - v)
}

deriv <- function(beta, a, lambda, penalty) {
  p <- length(beta)
  df <- numeric(p)
  if (penalty == "scad") {
    for (j in seq_len(p)) {
      if (abs(beta[j]) <= lambda) {
        df[j] <- lambda
      } else if (abs(beta[j]) <= a * lambda) {
        df[j] <- (a * lambda - abs(beta[j])) / (a - 1)
      }
    }
  } else if (penalty == "mcp") {
    for (j in seq_len(p)) {
      if (abs(beta[j]) <= a * lambda) df[j] <- lambda - abs(beta[j]) / a
    }
  } else {
    df <- rep(lambda, p)
  }
  df
}

QPADM <- function(X, y, tau, rho, lambda, iter = 500,
                  intercept = FALSE, penalty = "lasso", a = 3.7) {
  if (!(penalty %in% c("lasso", "scad", "mcp"))) penalty <- "lasso"
  n <- nrow(X)
  p <- ncol(X)
  if (intercept) {
    X <- cbind(1, X)
    beta <- rep(0, p + 1L)
  } else {
    beta <- rep(0, p)
  }
  u <- rep(0, n)
  r <- y - X %*% beta
  ABSTOL <- 1e-7
  RELTOL <- 1e-4
  iteration <- 1L
  while (iteration <= iter) {
    xbeta <- X %*% beta
    r <- shrinkcpp1(
      u / rho + y - xbeta - 0.5 * (2 * tau - 1) / (n * rho),
      0.5 / (n * rho)
    )
    beta_old <- beta
    df <- deriv(beta, a, lambda, penalty)
    for (j in seq_along(beta)) {
      xj <- X[, j]
      rj <- r + xj * beta[j]
      xjrj <- sum(xj * rj)
      xjsq <- sum(xj^2)
      beta[j] <- shrinkcpp1(xjrj / xjsq, df[j] / (rho * xjsq))
      r <- r + (beta_old[j] - beta[j]) * xj
    }
    u <- u + rho * (y - X %*% beta - r)
    rnorm <- sqrt(sum((y - X %*% beta - r)^2))
    snorm <- sqrt(sum((rho * X %*% (beta - beta_old))^2))
    epspri <- sqrt(n) * ABSTOL +
      RELTOL * max(c(sqrt(sum((X %*% beta)^2)), sqrt(sum(r^2))))
    epsdual <- sqrt(n) * ABSTOL + RELTOL * sqrt(sum(u^2))
    if (rnorm < epspri && snorm < epsdual) break
    iteration <- iteration + 1L
  }
  as.numeric(beta)
}

cv_qpadm_xtilde <- function(X, y, tau, lambda_grid, nfolds = 5,
                            rho = 1, intercept = FALSE,
                            penalty = "lasso", a = 3.7,
                            iter = 500, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  nfolds <- min(nfolds, n)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n))
  cv_mat <- matrix(NA_real_, nfolds, length(lambda_grid))
  for (k in seq_len(nfolds)) {
    id_val <- which(fold_id == k)
    id_tr <- setdiff(seq_len(n), id_val)
    Xtr <- X[id_tr, , drop = FALSE]
    ytr <- y[id_tr]
    Xva <- X[id_val, , drop = FALSE]
    yva <- y[id_val]
    for (j in seq_along(lambda_grid)) {
      b <- QPADM(
        Xtr, ytr, tau = tau, rho = rho, lambda = lambda_grid[j],
        iter = iter, intercept = intercept, penalty = penalty, a = a
      )
      cv_mat[k, j] <- check_loss(as.numeric(yva - Xva %*% b), tau)
    }
  }
  cv_mean <- colMeans(cv_mat)
  id_min <- which.min(cv_mean)
  list(lambda_hat = lambda_grid[id_min], cv_mean = cv_mean, id = id_min)
}

rqpen_get_beta <- function(rs, tau = 0.5, drop_intercept = TRUE) {
  if (!is.null(rs$coefficients)) {
    b <- as.numeric(rs$coefficients)
    if (drop_intercept && length(b) >= 2L) return(b[-1L])
    return(b)
  }
  if (!is.null(rs$models) && length(rs$models) > 0L) {
    nms <- names(rs$models)
    tau_pat <- paste0("tau", format(tau, scientific = FALSE, trim = TRUE))
    cand <- if (!is.null(nms) && length(nms) > 0L) {
      hit <- grep(tau_pat, nms, fixed = TRUE)
      if (length(hit) == 0L) seq_along(rs$models) else hit
    } else {
      seq_along(rs$models)
    }
    for (ii in cand) {
      obj <- rs$models[[ii]]
      if (!is.null(obj$coefficients)) {
        b <- as.numeric(obj$coefficients)
        if (drop_intercept && length(b) >= 2L) return(b[-1L])
        return(b)
      }
      if (!is.null(obj$model) && !is.null(obj$model$coefficients)) {
        b <- as.numeric(obj$model$coefficients)
        if (drop_intercept && length(b) >= 2L) return(b[-1L])
        return(b)
      }
    }
  }
  if (!is.null(rs$fit) && !is.null(rs$fit$coefficients)) {
    b <- as.numeric(rs$fit$coefficients)
    if (drop_intercept && length(b) >= 2L) return(b[-1L])
    return(b)
  }
  if (!is.null(rs$beta)) {
    b <- as.numeric(rs$beta)
    if (drop_intercept && length(b) >= 2L) return(b[-1L])
    return(b)
  }
  stop("Cannot locate coefficients in rqPen fit.", call. = FALSE)
}

cv_rqpen_xtilde <- function(X, y, lambda_grid, tau = 0.5,
                            nfolds = 5, alg = "br",
                            penalty = "LASSO", seed = NULL) {
  if (!has_pkg("rqPen")) {
    stop("Package rqPen is required for this competing method.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  nfolds <- min(nfolds, n)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n))
  cv_mat <- matrix(NA_real_, nfolds, length(lambda_grid))
  for (k in seq_len(nfolds)) {
    id_val <- which(fold_id == k)
    id_tr <- setdiff(seq_len(n), id_val)
    Xtr <- X[id_tr, , drop = FALSE]
    ytr <- y[id_tr]
    Xva <- X[id_val, , drop = FALSE]
    yva <- y[id_val]
    for (j in seq_along(lambda_grid)) {
      rs <- rqPen::rq.pen(y = as.vector(ytr), x = Xtr, alg = alg,
                          penalty = penalty, tau = tau,
                          lambda = lambda_grid[j])
      beta_hat <- rqpen_get_beta(rs, tau = tau)
      cv_mat[k, j] <- check_loss(as.numeric(yva - Xva %*% beta_hat), tau)
    }
  }
  cv_mean <- colMeans(cv_mat)
  cv_mean <- ifelse(is.na(cv_mean), 0, cv_mean)
  id_min <- which.min(cv_mean)
  list(lambda_hat = lambda_grid[id_min], id = id_min,
       cv_mean = cv_mean, cv_mat = cv_mat)
}

cv_qicd_xtilde <- function(X, y, lambda_grid, tau = 0.5,
                           nfolds = 5, intercept = FALSE,
                           funname = "lasso", seed = NULL) {
  if (!has_pkg("QICD")) {
    stop("Package QICD is required for this competing method.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  nfolds <- min(nfolds, n)
  fold_id <- sample(rep(seq_len(nfolds), length.out = n))
  cv_mat <- matrix(NA_real_, nfolds, length(lambda_grid))
  for (k in seq_len(nfolds)) {
    id_val <- which(fold_id == k)
    id_tr <- setdiff(seq_len(n), id_val)
    Xtr <- X[id_tr, , drop = FALSE]
    ytr <- y[id_tr]
    Xva <- X[id_val, , drop = FALSE]
    yva <- y[id_val]
    for (j in seq_along(lambda_grid)) {
      fit <- QICD::QICD(y = as.vector(ytr), x = Xtr, funname = funname,
                        tau = tau, lambda = lambda_grid[j],
                        intercept = intercept)
      beta_hat <- as.numeric(fit$beta_final)
      cv_mat[k, j] <- check_loss(as.numeric(yva - Xva %*% beta_hat), tau)
    }
  }
  cv_mean <- colMeans(cv_mat)
  cv_mean <- ifelse(is.na(cv_mean), 0, cv_mean)
  id_min <- which.min(cv_mean)
  list(lambda_hat = lambda_grid[id_min], id = id_min,
       cv_mean = cv_mean, cv_mat = cv_mat)
}

padmm_fun <- function() {
  if (exists("padmmR", mode = "function", inherits = TRUE)) {
    return(get("padmmR", mode = "function", inherits = TRUE))
  }
  if (has_pkg("openFHDQR") && "padmmR" %in% getNamespaceExports("openFHDQR")) {
    return(getExportedValue("openFHDQR", "padmmR"))
  }
  NULL
}

cv_padmm_xtilde <- function(X, y, tau, l1_seq, K = 5, seed = NULL) {
  padmmR <- padmm_fun()
  if (is.null(padmmR)) stop("padmmR is unavailable.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  K <- min(K, n)
  folds <- sample(rep(seq_len(K), length.out = n))
  val_loss <- numeric(length(l1_seq))
  for (i in seq_along(l1_seq)) {
    loss_i <- 0
    for (k in seq_len(K)) {
      train_idx <- which(folds != k)
      val_idx <- which(folds == k)
      Xt <- X[train_idx, , drop = FALSE]
      yt <- y[train_idx]
      Xv <- X[val_idx, , drop = FALSE]
      yv <- y[val_idx]
      fit <- padmmR(
        beta0 = rep(0, ncol(Xt)),
        z0 = yt - Xt %*% rep(0, ncol(Xt)),
        theta0 = rep(1, nrow(Xt)),
        sigma = 0.01,
        X = Xt,
        eta = eigen(crossprod(Xt), only.values = TRUE)$values[1],
        y = yt,
        l1 = l1_seq[i],
        l2 = 0,
        w = rep(1, ncol(Xt)),
        nu = rep(1, ncol(Xt)),
        tau = tau,
        gamma = 0.1
      )
      loss_i <- loss_i + check_loss(as.numeric(yv - Xv %*% fit$beta), tau)
    }
    val_loss[i] <- loss_i
  }
  best_idx <- which.min(val_loss)
  list(lambda.min = l1_seq[best_idx], val_loss = val_loss)
}

make_half_pairs <- function(X_raw, y_raw, permu) {
  n <- nrow(X_raw)
  p <- ncol(X_raw)
  stopifnot(length(permu) == n, n %% 2L == 0L)
  X_tilde <- matrix(0, nrow = n / 2L, ncol = p)
  y_tilde <- numeric(n / 2L)
  for (i in seq_len(n / 2L)) {
    i2 <- permu[2L * i]
    i1 <- permu[2L * i - 1L]
    X_tilde[i, ] <- X_raw[i2, ] - X_raw[i1, ]
    y_tilde[i] <- y_raw[i2] - y_raw[i1]
  }
  list(X = X_tilde, y = y_tilde)
}

subsample_l1_from_raw <- function(X_raw, y_raw, tau, lambda_grid,
                                  B = 10, nfolds = 5, seed = NULL) {
  if (!has_pkg("rqPen")) {
    stop("Package rqPen is required for this competing method.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X_raw)
  p <- ncol(X_raw)
  betas <- matrix(NA_real_, nrow = p, ncol = B)
  cvscores <- rep(NA_real_, B)
  lambdas <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    til <- make_half_pairs(X_raw, y_raw, sample.int(n))
    cvb <- cv_rqpen_xtilde(
      til$X, til$y, tau = tau, lambda_grid = lambda_grid,
      nfolds = nfolds, alg = "br", penalty = "LASSO",
      seed = if (is.null(seed)) NULL else seed + 1000L + b
    )
    rs <- rqPen::rq.pen(y = as.vector(til$y), x = til$X, alg = "br",
                        penalty = "LASSO", tau = tau,
                        lambda = cvb$lambda_hat)
    betas[, b] <- rqpen_get_beta(rs, tau = tau)
    cvscores[b] <- min(cvb$cv_mean)
    lambdas[b] <- cvb$lambda_hat
  }
  bbest <- which.min(cvscores)
  list(beta_hat = betas[, bbest], lambda_hat = lambdas[bbest],
       cvscore = cvscores[bbest], all_cv = cvscores,
       all_lambda = lambdas)
}

run_example1_competitors <- function(rep_id, seed, prep) {
  X <- prep$x
  y <- prep$y
  X_raw <- prep$x_original
  y_raw <- prep$y_original
  out <- list()

  out[["subsample_cv_xtilde"]] <- if (has_pkg("rqPen")) {
    safe_run({
      fit <- subsample_l1_from_raw(
        X_raw, y_raw, tau = 0.5, lambda_grid = seq(0.1, 0.3, by = 0.01),
        B = 10, nfolds = 5, seed = 15000 + rep_id
      )
      make_row(rep_id, seed, "subsample_cv_xtilde",
               fit$lambda_hat, fit$beta_hat)
    }, rep_id, seed, "subsample_cv_xtilde")
  } else {
    na_row(rep_id, seed, "subsample_cv_xtilde_missing")
  }

  out[["qpadm"]] <- safe_run({
    cvp <- cv_qpadm_xtilde(
      X, y, tau = 0.5, lambda_grid = seq(1000, 1500, by = 20),
      nfolds = 5, rho = 1, intercept = FALSE,
      penalty = "lasso", iter = 500, seed = 10000 + rep_id
    )
    b <- QPADM(X, y, tau = 0.5, rho = 1, lambda = cvp$lambda_hat,
               iter = 500, intercept = FALSE, penalty = "lasso")
    make_row(rep_id, seed, "qpadm_cv_xtilde", cvp$lambda_hat, b)
  }, rep_id, seed, "qpadm_cv_xtilde")

  out[["padmm"]] <- if (!is.null(padmm_fun())) {
    safe_run({
      padmmR <- padmm_fun()
      cvd <- cv_padmm_xtilde(
        X, y, tau = 0.5, l1_seq = lamb_list_resid, K = 5,
        seed = 12000 + rep_id
      )
      fitd <- padmmR(
        beta0 = rep(0, ncol(X)),
        z0 = y - X %*% rep(0, ncol(X)),
        theta0 = rep(1, nrow(X)),
        sigma = 0.01,
        X = X,
        eta = eigen(crossprod(X), only.values = TRUE)$values[1],
        y = y,
        l1 = cvd$lambda.min,
        l2 = 0,
        w = rep(1, ncol(X)),
        nu = rep(1, ncol(X)),
        tau = 0.5,
        gamma = 0.1
      )
      make_row(rep_id, seed, "padmm_cv_xtilde",
               cvd$lambda.min, as.numeric(fitd$beta))
    }, rep_id, seed, "padmm_cv_xtilde")
  } else {
    na_row(rep_id, seed, "padmm_cv_xtilde_missing")
  }

  out[["conquer"]] <- if (has_pkg("conquer")) {
    safe_run({
      fitc <- conquer::conquer.cv.reg(
        X, y, tau = 0.5, penalty = "lasso",
        lambdaSeq = seq(0.1, 0.6, by = 0.025)
      )
      make_row(rep_id, seed, "conquer_cv_xtilde",
               fitc$lambda.min, as.numeric(fitc$coeff.min[-1L]))
    }, rep_id, seed, "conquer_cv_xtilde")
  } else {
    na_row(rep_id, seed, "conquer_cv_xtilde_missing")
  }

  out[["hqreg"]] <- if (has_pkg("hqreg")) {
    safe_run({
      fith <- tryCatch(
        hqreg::cv.hqreg(X, y, method = "quantile", tau = 0.5,
                        lambda = seq(0.1, 0.6, by = 0.025), seed = rep_id),
        error = function(e) hqreg::cv.hqreg(X, y, method = "quantile",
                                            tau = 0.5, seed = rep_id)
      )
      btmp <- as.numeric(stats::coef(fith, lambda = "lambda.1se"))
      bh <- if (length(btmp) == ncol(X) + 1L) btmp[-1L] else btmp
      make_row(rep_id, seed, "hqreg_cv_xtilde", fith$lambda.1se, bh)
    }, rep_id, seed, "hqreg_cv_xtilde")
  } else {
    na_row(rep_id, seed, "hqreg_cv_xtilde_missing")
  }

  out[["qicd"]] <- if (has_pkg("QICD")) {
    safe_run({
      cvq <- cv_qicd_xtilde(
        X, y, lambda_grid = seq(250, 300, by = 5),
        tau = 0.5, nfolds = 5, intercept = FALSE,
        funname = "lasso", seed = 13000 + rep_id
      )
      fitq <- QICD::QICD(
        y = as.vector(y), x = X, funname = "lasso",
        tau = 0.5, lambda = cvq$lambda_hat, intercept = FALSE
      )
      make_row(rep_id, seed, "qicd_cv_xtilde",
               cvq$lambda_hat, as.numeric(fitq$beta_final))
    }, rep_id, seed, "qicd_cv_xtilde")
  } else {
    na_row(rep_id, seed, "qicd_cv_xtilde_missing")
  }

  out[["rqpen"]] <- if (has_pkg("rqPen")) {
    safe_run({
      cvr <- cv_rqpen_xtilde(
        X, y, lambda_grid = seq(0.1, 0.3, by = 0.01),
        tau = 0.5, nfolds = 5, alg = "br",
        penalty = "LASSO", seed = 14000 + rep_id
      )
      rs <- rqPen::rq.pen(
        y = as.vector(y), x = X, alg = "br",
        penalty = "LASSO", tau = 0.5, lambda = cvr$lambda_hat
      )
      make_row(rep_id, seed, "rqpen_cv_xtilde",
               cvr$lambda_hat, rqpen_get_beta(rs, tau = 0.5))
    }, rep_id, seed, "rqpen_cv_xtilde")
  } else {
    na_row(rep_id, seed, "rqpen_cv_xtilde_missing")
  }

  do.call(rbind, out)
}

run_one <- function(rep_id, seed) {
  set.seed(seed)
  dat <- simulate_rank_example(
    n = n0, p = p0, rho = rho0, error = "normal", beta = beta0
  )

  prep <- prepare_rank_data(dat$x, dat$y, standardize = TRUE,
                            rank_transform = TRUE)
  cold_start_beta <- glmnet_start(prep$x_original, prep$y_original)

  ## HBIC corresponds to the tuning step in the original Example 1 script.
  sel <- rank_admm_select_lambda(
    prep$x, prep$y,
    lambda_grid = lamb_list_resid,
    penalty = "L1",
    M = M0, G = G0, phi = phi0, rho = rho_admm,
    max_iter = HBIC_MAXITER_res,
    e.abs = HBIC_EABS, e.rel = HBIC_EABS,
    stop_rule = "residual",
    beta_init = cold_start_beta,
    rank_transform = FALSE,
    standardize = FALSE,
    n_original = n0,
    tau = 0.5, K = K0,
    freeze_phi = TRUE,
    residual_gated = TRUE,
    obj_scale = OBJ_SCALE,
    penalize_intercept = TRUE
  )

  ## Final reported estimator uses max_iter = 400 on the selected lambda.
  fit <- rank_admm_fit(
    prep$x, prep$y,
    penalty = "L1",
    lambda = sel$best_lambda,
    M = M0, G = G0, phi = phi0, rho = rho_admm,
    max_iter = FINAL_MAXITER_res,
    e.abs = HBIC_EABS, e.rel = HBIC_EABS,
    stop_rule = "residual",
    beta_start = cold_start_beta,
    rank_transform = FALSE,
    standardize = FALSE,
    tau = 0.5, K = K0,
    freeze_phi = TRUE,
    residual_gated = TRUE,
    obj_scale = OBJ_SCALE,
    penalize_intercept = TRUE
  )

  out <- list(
    ours = make_row(rep_id, seed, "ours_l1__select500_final400",
                    sel$best_lambda, fit$beta)
  )

  if (run_competitors) {
    out$competitors <- run_example1_competitors(
      rep_id = rep_id,
      seed = seed,
      prep = prep
    )
  }

  do.call(rbind, out)
}

example1_L1_results <- do.call(
  rbind,
  Map(run_one, rep_id = seq_along(seeds), seed = seeds)
)

print(example1_L1_results)
cat("Mean raw L1 error:",
    mean(example1_L1_results$err_l1, na.rm = TRUE), "\n")
