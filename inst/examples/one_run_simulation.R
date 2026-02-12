################################################################################
## Full simulation script
################################################################################

suppressPackageStartupMessages({
  library(MASS)
  library(foreach)
  library(doParallel)
  library(quantreg)
  library(rqPen)
  library(conquer)
  library(hqreg)
  library(QICD)
  library(openFHDQR)
})

## NOTE: must register a parallel backend, e.g.:
##   library(doParallel)
##   cl <- makeCluster(parallel::detectCores() - 1)
##   registerDoParallel(cl)

## ---------------------------
## settings
## ---------------------------
USE_TILDE <- TRUE

n0 <- 60
p0 <- 200
rho0 <- 0.5
sigma0 <- 1
tau0 <- 0.5

beta0 <- rep(0, p0)
beta0[1:4] <- 1

lamb_list_resid <- seq(0.3, 0.5, by = 0.01)
lamb_list_other <- seq(0.3, 0.5, by = 0.01)

M0 <- 2
G0 <- 4
K0 <- 1
phi0 <- 1/400
rho_admm <- 1.618

HBIC_MAXITER_res <- 200
HBIC_MAXITER_other <- 200
HBIC_EABS <- 1e-3
OBJ_SCALE <- "mean"

stop_grid <- list(
  list(name = "resid_final", stop_rule = "residual",
       e.abs = 1e-2, e.rel = 1e-2, freeze_phi_final = TRUE),
  list(name = "obj_final",   stop_rule = "objective",
       epsF = 1e-6, freeze_phi_final = TRUE),
  list(name = "beta_final",  stop_rule = "beta",
       epsBeta = 1e-4, freeze_phi_final = TRUE)
)

R <- 500

## ---------------------------
## Required functions checks
## ---------------------------
.must_exist <- function(fname) {
  if (!exists(fname, mode = "function")) {
    stop(sprintf("Missing required function '%s' in your environment.", fname), call. = FALSE)
  }
}
.must_exist("creating_x")
.must_exist("creating_y")
.must_exist("choose_lambda_hbic")
.must_exist("admm_both_stop")
.must_exist("threshold_beta")
.must_exist("calc_fp_fn")
.must_exist("HBIC")

## ---------------------------
##Data generator
## ---------------------------

AR <- function(n, p, rho) {
  idx <- 0:(p-1)
  Sigma <- rho^abs(outer(idx, idx, "-"))
  MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
}

prep_xy_tilde <- function(X_raw, y_raw) {
  X_raw <- scale(X_raw, center = TRUE, scale = TRUE)
  y_raw <- as.numeric(y_raw - mean(y_raw))
  if (USE_TILDE) {
    list(X = creating_x(X_raw), y = creating_y(y_raw), X_raw = X_raw, y_raw = y_raw)
  } else {
    list(X = X_raw, y = y_raw, X_raw = X_raw, y_raw = y_raw)
  }
}

## ---------------------------
## helpers
## ---------------------------
check_loss <- function(r, tau) sum(r * (tau - (r < 0)))

metrics_from_beta <- function(beta_hat, beta0, thr) {
  beta_thr <- threshold_beta(beta_hat, thr)
  sel_info <- calc_fp_fn(beta_thr, beta0, thr = thr)
  list(
    size  = length(sel_info$sel),
    err_l1 = sum(abs(beta_thr - beta0)),
    FP = sel_info$fp,
    FN = sel_info$fn,
    beta_thr = beta_thr
  )
}

make_row <- function(rep_id, method, lambda, beta_hat, beta0, thr) {
  met <- metrics_from_beta(beta_hat, beta0, thr = thr)
  data.frame(
    rep = rep_id,
    method = method,
    lambda = lambda,
    size = met$size,
    err_l1 = met$err_l1,
    FP = met$FP,
    FN = met$FN,
    V1 = met$beta_thr[1], V2 = met$beta_thr[2], V3 = met$beta_thr[3], V4 = met$beta_thr[4],
    stringsAsFactors = FALSE
  )
}

na_row <- function(rep_id, method) {
  data.frame(
    rep = rep_id, method = method, lambda = NA_real_,
    size = NA_integer_, err_l1 = NA_real_, FP = NA_integer_, FN = NA_integer_,
    V1 = NA_real_, V2 = NA_real_, V3 = NA_real_, V4 = NA_real_,
    stringsAsFactors = FALSE
  )
}

.safe_run <- function(expr, rep_id, method) {
  tryCatch(expr, error = function(e) {
    message(sprintf("[SKIP] %s failed at rep=%s: %s", method, rep_id, e$message))
    na_row(rep_id, method)
  })
}

################################################################################
## Baselines
################################################################################

## QPADM (Yu, Lin, Wang)
shrinkcpp1 <- function(u, v) {
  (1 + sign(u - v)) / 2 * (u - v) - (1 + sign(-u - v)) / 2 * (-u - v)
}
deriv <- function(beta, a, lambda, penalty) {
  p <- length(beta)
  df <- rep(0, p)
  if (penalty == "scad") {
    for (j in 1:p) {
      if (abs(beta[j]) <= lambda) {
        df[j] <- lambda
      } else if (abs(beta[j]) <= a * lambda) {
        df[j] <- (a * lambda - abs(beta[j])) / (a - 1)
      }
    }
  } else if (penalty == "mcp") {
    for (j in 1:p) {
      if (abs(beta[j]) <= a * lambda) {
        df[j] <- lambda - abs(beta[j]) / a
      }
    }
  } else {
    df <- rep(lambda, p)
  }
  df
}

QPADM <- function(X, y, tau, rho, lambda, iter = 500, intercept = FALSE,
                  M = 1, penalty = "lasso", a = 3.7) {
  if (!(penalty %in% c("lasso", "scad", "mcp"))) {
    warning(paste0('Penalty "', penalty, '" is not supported, Lasso estimation is returned'))
    penalty <- "lasso"
  }
  n <- nrow(X)
  p <- ncol(X)
  
  if (intercept) {
    X <- cbind(1, X)
    beta <- rep(0, p + 1)
  } else {
    beta <- rep(0, p)
  }
  
  u <- rep(0, n)
  r <- y - X %*% beta
  betaold <- beta
  iteration <- 1
  ABSTOL <- 1e-7
  RELTOL <- 1e-4
  
  while (iteration <= iter) {
    xbeta <- X %*% beta
    r <- shrinkcpp1(u / rho + y - xbeta - 0.5 * (2 * tau - 1) / (n * rho), 0.5 / (n * rho))
    
    betaold <- beta
    df <- deriv(beta, a, lambda, penalty)
    
    for (j in 1:length(beta)) {
      xj <- X[, j]
      rj <- r + xj * beta[j]
      xjrj <- sum(xj * rj)
      xjsq <- sum(xj^2)
      beta[j] <- shrinkcpp1(xjrj / xjsq, df[j] / (rho * xjsq))
      r <- r + (betaold[j] - beta[j]) * xj
    }
    
    u <- u + rho * (y - X %*% beta - r)
    
    rnorm <- sqrt(sum((y - X %*% beta - r)^2))
    snorm <- sqrt(sum((rho * X %*% (beta - betaold))^2))
    epspri <- sqrt(n) * ABSTOL + RELTOL * max(c(sqrt(sum((X %*% beta)^2)), sqrt(sum(r^2))))
    epsdual <- sqrt(n) * ABSTOL + RELTOL * sqrt(sum(u^2))
    if (rnorm < epspri && snorm < epsdual) break
    iteration <- iteration + 1
  }
  beta
}

cv_qpadm_xtilde <- function(X, y, tau, lambda_grid, nfolds = 5,
                            rho = 1, intercept = FALSE, penalty = "lasso", a = 3.7, iter = 500,
                            seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  fold_id <- sample(rep(1:nfolds, length.out = n))
  cv_mat <- matrix(NA_real_, nfolds, length(lambda_grid))
  for (k in 1:nfolds) {
    id_val <- which(fold_id == k)
    id_tr  <- setdiff(seq_len(n), id_val)
    Xtr <- X[id_tr, , drop = FALSE]; ytr <- y[id_tr]
    Xva <- X[id_val, , drop = FALSE]; yva <- y[id_val]
    for (j in seq_along(lambda_grid)) {
      b <- QPADM(Xtr, ytr, tau = tau, rho = rho, lambda = lambda_grid[j],
                 iter = iter, intercept = intercept, penalty = penalty, a = a)
      beta_hat <- as.numeric(b)
      r_va <- as.numeric(yva - Xva %*% beta_hat)
      cv_mat[k, j] <- check_loss(r_va, tau)
    }
  }
  cv_mean <- colMeans(cv_mat)
  jmin <- which.min(cv_mean)
  list(lambda_hat = lambda_grid[jmin], cv_mean = cv_mean, id = jmin)
}

## FDHQR
cv_padmm_xtilde <- function(X, y, tau, l1_seq, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  folds <- sample(rep(1:K, length.out = n))
  val_loss <- numeric(length(l1_seq))
  for (i in seq_along(l1_seq)) {
    loss_i <- 0
    for (k in 1:K) {
      train_idx <- which(folds != k)
      val_idx   <- which(folds == k)
      Xt_train <- X[train_idx, , drop = FALSE]
      yt_train <- y[train_idx]
      Xt_val   <- X[val_idx, , drop = FALSE]
      yt_val   <- y[val_idx]
      
      fit <- padmmR(
        beta0  = rep(0, ncol(Xt_train)),
        z0     = yt_train - Xt_train %*% rep(0, ncol(Xt_train)),
        theta0 = rep(1, nrow(Xt_train)),
        sigma  = 0.01,
        X      = Xt_train,
        eta    = eigen(t(Xt_train) %*% Xt_train, only.values = TRUE)$values[1],
        y      = yt_train,
        l1     = l1_seq[i], l2 = 0,
        w      = rep(1, ncol(Xt_train)),
        nu     = rep(1, ncol(Xt_train)),
        tau    = tau,
        gamma  = 0.1
      )
      preds <- Xt_val %*% fit$beta
      loss_i <- loss_i + check_loss(yt_val - preds, tau)
    }
    val_loss[i] <- loss_i
  }
  best_idx <- which.min(val_loss)
  list(lambda.min = l1_seq[best_idx], val_loss = val_loss)
}

## subsample
make_tilde_once <- function(X_raw, y_raw, permu) {
  n <- nrow(X_raw); p <- ncol(X_raw)
  stopifnot(length(permu) == n)
  stopifnot(n %% 2 == 0)
  X_tilde <- matrix(0, nrow = n/2, ncol = p)
  y_tilde <- numeric(n/2)
  for (i in 1:(n/2)) {
    i2 <- permu[2*i]; i1 <- permu[2*i - 1]
    X_tilde[i, ] <- X_raw[i2, ] - X_raw[i1, ]
    y_tilde[i]   <- y_raw[i2]   - y_raw[i1]
  }
  list(X = X_tilde, y = y_tilde)
}

rqpen_beta <- function(rs, tau = 0.5) {
  # robust extraction
  if (!is.null(rs$models) && length(rs$models) >= 1 &&
      !is.null(rs$models[[1]]$coefficients)) {
    return(as.numeric(rs$models[[1]]$coefficients[-1]))
  }
  if (!is.null(rs$coefficients)) return(as.numeric(rs$coefficients[-1]))
  stop("Cannot extract coefficients from rq.pen fit.")
}

cv_rqpen_on_tilde <- function(X_tilde, y_tilde, tau, lambda_grid, nfolds = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X_tilde)
  fold_id <- sample(rep(1:nfolds, length.out = n))
  cv_mat <- matrix(NA_real_, nfolds, length(lambda_grid))
  
  for (k in 1:nfolds) {
    id_val <- which(fold_id == k)
    id_tr  <- setdiff(seq_len(n), id_val)
    Xtr <- X_tilde[id_tr, , drop = FALSE]; ytr <- y_tilde[id_tr]
    Xva <- X_tilde[id_val, , drop = FALSE]; yva <- y_tilde[id_val]
    
    for (j in seq_along(lambda_grid)) {
      rs <- rq.pen(y = as.vector(ytr), x = Xtr, alg = "br",
                   penalty = "LASSO", tau = tau, lambda = lambda_grid[j])
      b <- rqpen_beta(rs, tau)
      rva <- as.numeric(yva - Xva %*% b)
      cv_mat[k, j] <- check_loss(rva, tau)
    }
  }
  
  cv_mean <- colMeans(cv_mat)
  jbest <- which.min(cv_mean)
  list(lambda_hat = lambda_grid[jbest], cv_mean = cv_mean, id = jbest)
}

baseline_subsample_cv_xtilde <- function(X_raw, y_raw, tau = 0.5,
                                         lambda_grid,
                                         B = 10,
                                         nfolds = 5,
                                         agg = c("best_cv", "median"),
                                         seed = NULL) {
  agg <- match.arg(agg)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X_raw); p <- ncol(X_raw)
  stopifnot(n %% 2 == 0)
  
  betas <- matrix(NA_real_, nrow = p, ncol = B)
  cvscores <- rep(NA_real_, B)
  lambdas  <- rep(NA_real_, B)
  
  for (b in 1:B) {
    permu <- sample.int(n)
    til <- make_tilde_once(X_raw, y_raw, permu)
    X_t <- til$X; y_t <- til$y
    
    cvb <- cv_rqpen_on_tilde(X_t, y_t, tau = tau, lambda_grid = lambda_grid,
                             nfolds = nfolds, seed = if (is.null(seed)) NULL else seed + 1000 + b)
    
    rs <- rq.pen(y = as.vector(y_t), x = X_t, alg = "br",
                 penalty = "LASSO", tau = tau, lambda = cvb$lambda_hat)
    betas[, b] <- rqpen_beta(rs, tau)
    cvscores[b] <- min(cvb$cv_mean)
    lambdas[b]  <- cvb$lambda_hat
  }
  
  if (agg == "best_cv") {
    bbest <- which.min(cvscores)
    list(beta_hat = betas[, bbest], lambda_hat = lambdas[bbest],
         cvscore = cvscores[bbest], all_cv = cvscores, all_lambda = lambdas)
  } else {
    list(beta_hat = apply(betas, 1, median), lambda_hat = median(lambdas),
         cvscore = median(cvscores), all_cv = cvscores, all_lambda = lambdas)
  }
}



## QICD
cv_qicd_xtilde <- function(X, y, lambda_grid,
                           tau = 0.5, nfolds = 5, intercept = FALSE,
                           funname = "lasso", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  fold_id <- sample(rep(1:nfolds, length.out = n))
  cv_mat <- matrix(NA_real_, nfolds, length(lambda_grid))
  for (k in 1:nfolds) {
    id_val <- which(fold_id == k)
    id_tr  <- setdiff(seq_len(n), id_val)
    X_tr <- X[id_tr, , drop = FALSE]; y_tr <- y[id_tr]
    X_va <- X[id_val, , drop = FALSE]; y_va <- y[id_val]
    for (j in seq_along(lambda_grid)) {
      fit <- QICD(y = as.vector(y_tr), x = X_tr, funname = funname,
                  tau = tau, lambda = lambda_grid[j], intercept = intercept)
      beta_hat <- as.numeric(fit$beta_final)
      r_va <- as.numeric(y_va - X_va %*% beta_hat)
      cv_mat[k, j] <- check_loss(r_va, tau)
    }
  }
  cv_mean <- colMeans(cv_mat)
  id_min <- which.min(cv_mean)
  list(lambda_hat = lambda_grid[id_min], id = id_min, cv_mean = cv_mean, cv_mat = cv_mat)
}

## ---- rqPen CV ----
rqpen_get_beta <- function(rs, tau = 0.5) {
  key <- paste0("tau", tau, "a1")
  if (!is.null(rs$models) && !is.null(rs$models[[key]]) &&
      !is.null(rs$models[[key]]$coefficients)) {
    b <- rs$models[[key]]$coefficients
    return(as.numeric(b[-1]))
  }
  if (!is.null(rs$coefficients)) {
    b <- rs$coefficients
    return(as.numeric(b[-1]))
  }
  stop("Cannot find coefficients in rq.pen object. Inspect rs$models names.")
}

cv_rqpen_xtilde <- function(X, y, lambda_grid,
                            tau = 0.5, nfolds = 5, alg = "br",
                            penalty = "LASSO", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  fold_id <- sample(rep(1:nfolds, length.out = n))
  cv_mat <- matrix(NA_real_, nfolds, length(lambda_grid))
  for (k in 1:nfolds) {
    id_val <- which(fold_id == k)
    id_tr  <- setdiff(seq_len(n), id_val)
    X_tr <- X[id_tr, , drop = FALSE]; y_tr <- y[id_tr]
    X_va <- X[id_val, , drop = FALSE]; y_va <- y[id_val]
    for (j in seq_along(lambda_grid)) {
      rs <- rq.pen(y = as.vector(y_tr), x = X_tr, alg = alg, penalty = penalty,
                   tau = tau, lambda = lambda_grid[j])
      beta_hat <- rqpen_get_beta(rs, tau = tau)
      r_va <- as.numeric(y_va - X_va %*% beta_hat)
      cv_mat[k, j] <- check_loss(r_va, tau)
    }
  }
  cv_mean <- colMeans(cv_mat)
  id_min <- which.min(cv_mean)
  list(lambda_hat = lambda_grid[id_min], id = id_min, cv_mean = cv_mean, cv_mat = cv_mat)
}

run_others_one_rep <- function(r, X_raw, y_raw, beta0, tau0,
                               lamb_list_qpadm = lamb_list_other,
                               lamb_list_padmm = lamb_list_other,
                               lamb_list_qicd  = lamb_list_other,
                               lamb_list_rqpen = lamb_list_other,
                               lamb_list_sub   = lamb_list_other,
                               J_sub = 20,
                               frac_sub = 0.5,
                               lamb_list_conquer = NULL,
                               lamb_list_hqreg = NULL,
                               seed_base = 10000,
                               BETA_THR_OTHERS = 0) {
  xy <- prep_xy_tilde(X_raw, y_raw)
  X <- xy$X; y <- xy$y
  
  out <- list()
  
  out[["subsample_cv_xtilde"]] <- .safe_run({
    Xraw_std <- xy$X_raw
    yraw_ctr <- xy$y_raw
    
    fit <- baseline_subsample_cv_xtilde(
      X_raw = Xraw_std, y_raw = yraw_ctr,
      tau = tau0,
      lambda_grid = lamb_list_rqpen, 
      B = 2*5, 
      nfolds = 5,
      agg = "best_cv",
      seed = seed_base + 5000 + r
    )
    make_row(r, "subsample_cv_xtilde", fit$lambda_hat,
             beta_hat = fit$beta_hat, beta0 = beta0, thr = BETA_THR_OTHERS)
  }, r, "subsample_cv_xtilde")
  
  
  # QPADM
  out[["qpadm"]] <- .safe_run({
    cvp <- cv_qpadm_xtilde(X, y, tau = tau0, lambda_grid = lamb_list_qpadm,
                           nfolds = 5, rho = 1, intercept = FALSE,
                           penalty = "lasso", iter = 500, seed = seed_base + r)
    b <- QPADM(X, y, tau = tau0, rho = 1, lambda = cvp$lambda_hat,
               iter = 500, intercept = FALSE, penalty = "lasso")
    make_row(r, "qpadm_cv_xtilde", cvp$lambda_hat, beta_hat = as.numeric(b), beta0 = beta0, thr = BETA_THR_OTHERS)
  }, r, "qpadm_cv_xtilde")
  
  # padmmR
  out[["padmm"]] <- if (exists("padmmR", mode = "function")) .safe_run({
    cvd <- cv_padmm_xtilde(X, y, tau = tau0, l1_seq = lamb_list_padmm, K = 5, seed = seed_base + 2000 + r)
    fitd <- padmmR(
      beta0  = rep(0, ncol(X)),
      z0     = y - X %*% rep(0, ncol(X)),
      theta0 = rep(1, nrow(X)),
      sigma  = 0.01,
      X      = X,
      eta    = eigen(t(X) %*% X, only.values = TRUE)$values[1],
      y      = y,
      l1     = cvd$lambda.min, l2 = 0,
      w      = rep(1, ncol(X)),
      nu     = rep(1, ncol(X)),
      tau    = tau0,
      gamma  = 0.1
    )
    make_row(r, "padmm_cv_xtilde", cvd$lambda.min, beta_hat = as.numeric(fitd$beta), beta0 = beta0, thr = BETA_THR_OTHERS)
  }, r, "padmm_cv_xtilde") else na_row(r, "padmm_cv_xtilde")
  
  # conquer
  out[["conquer"]] <- if (exists("conquer.cv.reg", mode = "function")) .safe_run({
    fitc <- conquer.cv.reg(X, y, tau = tau0, penalty = "lasso", lambdaSeq = lamb_list_conquer)
    make_row(r, "conquer_cv_xtilde", fitc$lambda.min, beta_hat = as.numeric(fitc$coeff.min[-1]), beta0 = beta0, thr = BETA_THR_OTHERS)
  }, r, "conquer_cv_xtilde") else na_row(r, "conquer_cv_xtilde")
  
  # hqreg
  out[["hqreg"]] <- if (exists("cv.hqreg", mode = "function")) .safe_run({
    
    
    fith <- tryCatch(
      cv.hqreg(X, y, method = "quantile", tau = tau0, lambda = lamb_list_hqreg, seed = r),
      error = function(e) cv.hqreg(X, y, method = "quantile", tau = tau0, seed = r)  # fallback to internal grid
    )
    
    btmp <- as.numeric(coef(fith, lambda = "lambda.1se"))
    
    # drop intercept only if present
    if (length(btmp) == ncol(X) + 1) bh <- btmp[-1] else bh <- btmp
    
    make_row(r, "hqreg_cv_xtilde", fith$lambda.1se, beta_hat = bh, beta0 = beta0, thr = BETA_THR_OTHERS)
    
  }, r, "hqreg_cv_xtilde") else na_row(r, "hqreg_cv_xtilde")
  
  # QICD
  out[["qicd"]] <- if (exists("QICD", mode = "function")) .safe_run({
    cvq <- cv_qicd_xtilde(X, y, lambda_grid = lamb_list_qicd,
                          tau = tau0, nfolds = 5, intercept = FALSE,
                          funname = "lasso", seed = seed_base + 3000 + r)
    fitq <- QICD(y = as.vector(y), x = X, funname = "lasso", tau = tau0,
                 lambda = cvq$lambda_hat, intercept = FALSE)
    make_row(r, "qicd_cv_xtilde", cvq$lambda_hat, beta_hat = as.numeric(fitq$beta_final), beta0 = beta0, thr = BETA_THR_OTHERS)
  }, r, "qicd_cv_xtilde") else na_row(r, "qicd_cv_xtilde")
  
  # rqPen
  out[["rqpen"]] <- if (exists("rq.pen", mode = "function")) .safe_run({
    cvr <- cv_rqpen_xtilde(X, y, lambda_grid = lamb_list_rqpen,
                           tau = tau0, nfolds = 5, alg = "br",
                           penalty = "LASSO", seed = seed_base + 4000 + r)
    rs <- rq.pen(y = as.vector(y), x = X, alg = "br", penalty = "LASSO",
                 tau = tau0, lambda = cvr$lambda_hat)
    br <- rqpen_get_beta(rs, tau = tau0)
    make_row(r, "rqpen_cv_xtilde", cvr$lambda_hat, beta_hat = br, beta0 = beta0, thr = BETA_THR_OTHERS)
  }, r, "rqpen_cv_xtilde") else na_row(r, "rqpen_cv_xtilde")
  
  do.call(rbind, out)
}


run_pipeline_all <- function() {
  
  seed_master <- 2021
  dat_list <- lapply(1:R, function(r) {
    set.seed(seed_master + r)
    X_raw <- AR(n0, p0, rho0)
    y_raw <- as.numeric(X_raw %*% beta0 + rnorm(n0, 0, sigma0))
    list(X_raw = X_raw, y_raw = y_raw)
  })
  
  
  res_ours <- foreach(
    r = 1:R,
    .combine = rbind,
    .packages = c("MASS", "foreach", "quantreg", "glmnet"),
    .export = c("AR","creating_x","creating_y","matsplit_mg","prox_beta",
                "l2norm","support_set","calc_fp_fn",
                "threshold_beta","admm_both_stop","HBIC","choose_lambda_hbic",
                "prep_xy_tilde",
                "n0","p0","rho0","sigma0","tau0","beta0",
                "USE_TILDE","lamb_list_resid","lamb_list_other",
                "M0","G0","K0","phi0","rho_admm",
                "HBIC_MAXITER_res","HBIC_MAXITER_other","HBIC_EABS","OBJ_SCALE",
                "stop_grid")
  ) %dopar% {
    
    seed_master <- 2021
    set.seed(seed_master + as.integer(r))
    
    X_raw <- dat_list[[r]]$X_raw
    y_raw <- dat_list[[r]]$y_raw
    
    
    X_raw <- scale(X_raw, center = TRUE, scale = TRUE)
    y_raw <- as.numeric(y_raw - mean(y_raw))
    
    if (USE_TILDE) {
      X <- creating_x(X_raw)
      y <- creating_y(y_raw)
      x_ori_for_pen <- X_raw
    } else {
      X <- X_raw
      y <- y_raw
      x_ori_for_pen <- X_raw
    }
    
    # cold start
    cv_lasso <- cv.glmnet(X_raw, y_raw, alpha = 1, intercept = TRUE)
    beta_lasso_raw <- as.numeric(coef(cv_lasso, s = "lambda.min"))
    cold_start_beta <- beta_lasso_raw[-1]
    
    # HBIC selection
    sel_resid <- choose_lambda_hbic(
      lamb_list = lamb_list_resid, X = X, y = y, x.ori = x_ori_for_pen, beta_init_guess = cold_start_beta,
      stop_rule = "residual", freeze_phi = TRUE, tau = tau0, M = M0, G = G0, K = K0,
      phi = phi0, rho = rho_admm, max_iter = HBIC_MAXITER_res, e.abs = HBIC_EABS, e.rel = HBIC_EABS,
      residual_gated = TRUE, obj_scale = OBJ_SCALE
    )
    
    sel_other <- choose_lambda_hbic(
      lamb_list = lamb_list_other, X = X, y = y, x.ori = x_ori_for_pen, beta_init_guess = cold_start_beta,
      stop_rule = "objective", freeze_phi = TRUE, tau = tau0, M = M0, G = G0, K = K0,
      phi = phi0, rho = rho_admm, max_iter = HBIC_MAXITER_other, e.abs = HBIC_EABS, e.rel = HBIC_EABS,
      residual_gated = TRUE, obj_scale = OBJ_SCALE
    )
    
    out_rows <- vector("list", length(stop_grid))
    
    for (j in seq_along(stop_grid)) {
      cfg <- stop_grid[[j]]
      
      if (cfg$stop_rule == "residual") {
        current_lambda <- sel_resid$best_lambda
        current_start  <- sel_resid$best_model$beta.hat
      } else {
        current_lambda <- sel_other$best_lambda
        current_start  <- sel_other$best_model$beta.hat
      }
      
      eabs_j   <- if (!is.null(cfg$e.abs))   cfg$e.abs   else HBIC_EABS
      erel_j   <- if (!is.null(cfg$e.rel))   cfg$e.rel   else HBIC_EABS
      epsF_j   <- if (!is.null(cfg$epsF))    cfg$epsF    else 1e-6
      epsB_j   <- if (!is.null(cfg$epsBeta)) cfg$epsBeta else 1e-4
      freeze_j <- if (!is.null(cfg$freeze_phi_final)) cfg$freeze_phi_final else TRUE
      
      fit <- admm_both_stop(
        X = X, y = y,
        lambda = current_lambda,
        beta_start = current_start,
        max_iter = 500,
        rho = rho_admm,
        tau = tau0,
        e.abs = eabs_j,
        e.rel = erel_j,
        phi = phi0,
        M = M0,
        G = G0,
        K = K0,
        stop_rule = cfg$stop_rule,
        epsF = epsF_j,
        epsBeta = epsB_j,
        delta_support = 1e-6,
        L_support = 5,
        penalize_intercept = TRUE,
        min_iter = 10,
        residual_gated = TRUE,
        obj_scale = OBJ_SCALE,
        freeze_phi = freeze_j,
        verbose = FALSE
      )
      
      
      BETA_THR <- 0.2
      BETA_THR_OTHERS <- 0
      beta_hat_thr <- threshold_beta(fit$beta.hat, BETA_THR)
      err_l1 <- sum(abs(beta_hat_thr - beta0))
      sel_info <- calc_fp_fn(beta_hat_thr, beta0, thr = BETA_THR)
      
      metrics_row <- data.frame(
        rep = r,
        method = paste0("ours__", cfg$name),
        lambda = current_lambda,
        size = length(sel_info$sel),
        err_l1 = err_l1,
        FP = sel_info$fp,
        FN = sel_info$fn,
        stringsAsFactors = FALSE
      )
      
      beta_row <- as.data.frame(t(beta_hat_thr[1:4]))
      colnames(beta_row) <- paste0("V", 1:4)
      
      out_rows[[j]] <- cbind(metrics_row, beta_row)
    }
    
    do.call(rbind, out_rows)
  }
  
  # baselines serial (all xy_tilde)
  
  res_others <- do.call(rbind, lapply(1:R, function(r) {
    dat <- dat_list[[r]]
    run_others_one_rep(
      r = r,
      X_raw = dat$X_raw,
      y_raw = dat$y_raw,
      beta0 = beta0,
      tau0 = tau0,
      lamb_list_qpadm = seq(1000, 1500, by=20),
      lamb_list_padmm = lamb_list_other,
      lamb_list_qicd  = seq(250, 300, by=5),
      lamb_list_rqpen = seq(0.1, 0.3, by=0.01),
      lamb_list_sub   = lamb_list_other,
      lamb_list_conquer = seq(0.1, 0.6, by = 0.025),
      lamb_list_hqreg = seq(0.1, 0.6, by = 0.025),
      J_sub = 10,
      frac_sub = 0.5,
      BETA_THR_OTHERS = 0
    )
  }))
  
  rbind(res_ours, res_others)
}

## ---------------------------
## 7) RUN
## ---------------------------
res_df <- run_pipeline_all()

## If you started a cluster manually:
## stopCluster(cl)
