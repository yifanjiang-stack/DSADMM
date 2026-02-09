USE_TILDE <- TRUE

n0 <- 60
p0 <- 200
rho0 <- 0.5
sigma0 <- 1
tau0 <- 0.5

beta0 <- rep(0, p0)
beta0[1:4] <- 1

# 1. LAMBDA SETTINGS
lamb_list_resid <- seq(0.01, 1, by = 0.01)
lamb_list_other <- seq(0.01, 1, by = 0.01)

M0 <- 1
G0 <- 8
K0 <- 1
phi0 <- 1/400
rho_admm <- 1.618

# 2. ITERATION LIMITS (Used during selection only now)
HBIC_MAXITER_res <- 200
HBIC_MAXITER_other <- 200
HBIC_EABS <- 1e-3
OBJ_SCALE <- "mean"
BETA_THR <- 0

# UPDATED: per-stop-rule settings; each row will be actually run under its own stop_rule
stop_grid <- list(
  list(name = "resid_final", stop_rule = "residual",
       e.abs = 1e-2, e.rel = 1e-2, freeze_phi_final = TRUE),
  list(name = "obj_final",   stop_rule = "objective",
       epsF = 1e-6, freeze_phi_final = TRUE),
  list(name = "beta_final",  stop_rule = "beta",
       epsBeta = 1e-4, freeze_phi_final = TRUE)
)

R <- 500

res_df <- foreach(
  r = 1:R,
  .combine = rbind,
  .packages = c("MASS", "foreach", "quantreg"),
  .export = c("AR","creating_x","creating_y","matsplit_mg","prox_beta",
              "l2norm","support_set","calc_fp_fn",
              "threshold_beta","admm_both_stop","HBIC","choose_lambda_hbic")
) %dopar% {

  X_raw <- AR(n0, p0, rho0)
  y_raw <- as.numeric(X_raw %*% beta0 + rnorm(n0, 0, sigma0))

  X_raw <- scale(X_raw, center = TRUE, scale = TRUE)
  y_raw <- y_raw - mean(y_raw)

  if (USE_TILDE) {
    X <- creating_x(X_raw)
    y <- creating_y(y_raw)
    x_ori_for_pen <- X_raw
  } else {
    X <- X_raw
    y <- y_raw
    x_ori_for_pen <- X_raw
  }

  library(glmnet)
  cv_lasso <- cv.glmnet(X_raw, y_raw, alpha = 1, intercept = TRUE)
  beta_lasso_raw <- as.numeric(coef(cv_lasso, s = "lambda.min"))
  cold_start_beta <- beta_lasso_raw[-1]

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

    beta_hat_thr <- threshold_beta(fit$beta.hat, BETA_THR)
    err_l1 <- sum(abs(beta_hat_thr - beta0))
    sel_info <- calc_fp_fn(beta_hat_thr, beta0, thr = 0)

    metrics_row <- data.frame(
      rep = r,
      stop_name = cfg$name,
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
