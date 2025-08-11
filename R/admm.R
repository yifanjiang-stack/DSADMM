#' Double-splitting ADMM for quantile regression
#'
#' @param X Design matrix (n x p)
#' @param y Response vector (length n)
#' @param lambda L1 penalty parameter
#' @param max_iter Maximum ADMM iterations
#' @param rho,phi ADMM parameters
#' @param tau Quantile level in (0,1)
#' @param e.rel,e.abs Relative/absolute tolerances
#' @param M,G,K Splitting/acceleration settings
#' @param Args Optional list of warm-starts: beta, z, w, eta
#' @param beta_true Optional true beta for diagnostics
#' @param idx_true Optional support set of true beta (indices of nonzeros)
#'
#' @return A list with elements
#' \itemize{
#'   \item beta.hat: estimated coefficient vector (length p)
#'   \item iter: number of iterations used
#'   \item err, FP, FN, size, P1, P2: diagnostics (NA if not provided true values)
#' }
#'
#' @import MASS
#' @import doParallel
#' @import foreach
#' @importFrom foreach %dopar% %do% foreach
#' @importFrom MASS ginv
#' @export
admm_both <- function(X, y, lambda = 0, max_iter, rho = 1, tau,
                      e.rel = 1e-3, e.abs = 1e-3, phi, M = 2, G = 2, K = 1,
                      Args = NULL, beta_true = NULL, idx_true = NULL){

  p <- ncol(X); n <- nrow(X)

  # Matrix partitions
  x.mg <- matsplitter(X, M, G)
  y.m  <- matrix(y, ncol = M)

  # Initialization
  if (is.null(Args)){
    z.m <- matrix(0, nrow = n/M, ncol = M)
    beta <- matrix(0, nrow = p, ncol = M)
    beta.mg <- matsplitter(beta, G, 1)
    beta.g <- matrix(0, nrow = p/G, ncol = G)
    w <- matrix(0, nrow = n, ncol = G)
    w.mg <- matsplitter(w, M, 1)
    gamma <- matrix(0, nrow = p, ncol = M)
    gamma.g <- matsplitter(gamma, G, 1)
    eta <- matrix(0, nrow = n, ncol = G)
    eta.mg <- matsplitter(eta, M, 1)
    theta.m <- matrix(0, nrow = n/M, ncol = M)
  } else{
    beta <- Args$beta
    z.m  <- Args$z
    w    <- Args$w
    eta  <- Args$eta
    beta.mg <- matsplitter(beta, G, 1)
    beta.g  <- matrix(0, nrow = p/G, ncol = G)
    w.mg    <- matsplitter(w, M, 1)
    gamma   <- matrix(0, nrow = p, ncol = M)
    gamma.g <- matsplitter(gamma, G, 1)
    eta.mg  <- matsplitter(eta, M, 1)
    theta.m <- matrix(0, nrow = n/M, ncol = M)
  }

  r_1 <- matsplitter(matrix(0, nrow = p, ncol = M), G, 1)
  r_2 <- matsplitter(matrix(0, nrow = n, ncol = G), M, 1)
  r_3 <- matrix(0, nrow = n/M, ncol = M)
  abeta <- matsplitter(matrix(0, nrow = n, ncol = G), M, 1)

  # Inverse step (parallelized)
  doParallel::registerDoParallel(M*G)
  x.inv <- foreach::foreach (i=1:(M*G)) %dopar% {
    MASS::ginv(t(x.mg[[i]]) %*% x.mg[[i]] + diag(ncol(x.mg[[i]])))
  }
  x.invnew <- matrix(x.inv, ncol = G)

  # Residuals
  r.norm <- 1; s.norm <- 1; e.primal <- 0; e.dual <- 0; t <- 1L

  # ADMM loop
  while ((r.norm > e.primal || s.norm > e.dual) && t < max_iter){

    z.1 <- z.m; w.1 <- w.mg; eta.1 <- eta.mg

    # beta.g
    doParallel::registerDoParallel(G)
    beta.g[,1:G] <- foreach::foreach(g=1:G, .combine=cbind) %dopar% {
      u_vec <- rowMeans(beta.mg[[g]]) - (1/phi)*rowMeans(gamma.g[[g]])
      a_vec <- (lambda)/(phi*M)
      prox_beta(u_vec, a_vec)
    }

    # beta.mg
    doParallel::registerDoParallel(G)
    beta.mg <- foreach::foreach(g=1:G) %dopar% {
      foreach::foreach(m=1:M, .combine=cbind) %do% {
        x.invnew[m,g][[1]] %*%
          (beta.g[,g] - 1/phi * t(x.mg[m,g][[1]]) %*% eta.mg[[m]][,g] +
             t(x.mg[m,g][[1]]) %*% w.mg[[m]][,g] + 1/phi * gamma.g[[g]][,m])
      }
    }

    # w.mg intermediate
    doParallel::registerDoParallel(G)
    for (m in 1:M){
      for (g in 1:G){
        abeta[[m]][,g] <- x.mg[m,g][[1]] %*% beta.mg[[g]][,m]
      }
    }

    w.mg <- foreach::foreach(m=1:M) %do% {
      foreach::foreach(g=1:G, .combine=cbind) %dopar% {
        1/G * (y.m[,m] - z.m[,m] + G * x.mg[m,g][[1]] %*% beta.mg[[g]][,m] -
                 rowSums(abeta[[m]]))
      }
    }

    # z.m
    doParallel::registerDoParallel(M)
    z.m[,1:M] <- foreach::foreach(m=1:M, .combine=cbind) %dopar% {
      xi_vec <- y.m[,m] - rowSums(w.mg[[m]]) - (theta.m[,m]/phi)
      alpha_vec <- 1/(n*phi/M)
      prox_z(xi = xi_vec, alpha = alpha_vec, quan = tau)
    }

    # duals
    for (g in 1:G){
      for (m in 1:M){
        gamma.g[[g]][,m] <- gamma.g[[g]][,m] + rho*phi*(-beta.mg[[g]][,m]+beta.g[,g])
      }
    }
    for (m in 1:M){
      theta.m[,m] <- theta.m[,m] + rho*phi*(z.m[,m] + rowSums(w.mg[[m]]) - y.m[,m])
    }
    for (m in 1:M){
      for (g in 1:G){
        eta.mg[[m]][,g] <- eta.mg[[m]][,g] +
          rho*phi*(-w.mg[[m]][,g] + x.mg[m,g][[1]] %*% beta.mg[[g]][,m])
      }
    }

    t <- t + 1L

    # residuals
    for (g in 1:G){
      for (m in 1:M){
        r_1[[g]][,m] <- -beta.mg[[g]][,m] + beta.g[,g]
      }
    }
    for (m in 1:M){
      for (g in 1:G){
        r_2[[m]][,g] <- -w.mg[[m]][,g] + x.mg[m,g][[1]] %*% beta.mg[[g]][,m]
      }
    }
    for (m in 1:M){
      r_3[,m] <- -y.m[,m] + z.m[,m] + rowSums(w.mg[[m]])
    }

    s <- phi/M * t(X) %*% (c(z.m) - c(z.1))

    r.norm <- l2norm(unlist(r_1)) + l2norm(unlist(r_2)) + l2norm(c(r_3))
    s.norm <- l2norm(s)

    e.primal <- sqrt(n) * e.abs +
      e.rel * max(l2norm(X%*%c(beta.g)), l2norm(c(z.m)), l2norm(y)) / sqrt(M*G)
    e.dual <- sqrt(p) * e.abs + e.rel * l2norm(t(X) %*% c(theta.m)) / (M)

    # adapt phi
    if ((t %% K) == 0){
      if (r.norm/M > 10*s.norm) {
        phi <- 2*phi
      } else if (s.norm > 10*r.norm/M) {
        phi <- phi/2
      }
    }
  }

  beta <- c(beta.g)

  # diagnostics (optional, only if true info provided)
  P1 <- P2 <- FP <- FN <- err <- NA_real_
  if (!is.null(idx_true)) {
    P1 <- as.numeric(sum(beta[idx_true] != 0) == length(idx_true))
    P2 <- as.numeric(beta[1] != 0)
    FP <- sum(beta[-idx_true] != 0) / length(beta[-idx_true])
    FN <- sum(beta[idx_true] == 0) / length(idx_true)
  }
  if (!is.null(beta_true)) {
    err <- l1norm(beta - beta_true)
  }

  list(beta.hat = beta,
       iter = t,
       err = err,
       P1 = P1,
       P2 = P2,
       FP = FP,
       FN = FN,
       size = sum(beta != 0))
}

#' HBIC-style score for model selection
#' @param lamb penalty value
#' @param X,y design and response
#' @param K,M,G,phi,rho,tau,e.rel,e.abs,max_iter algorithm settings (subset used)
#' @param x.ori original (un-transformed) X for n computation
#' @return list(hbic=..., model=...)
#' @export
HBIC <- function(lamb, X, y, K = 1, x.ori, M = 2, G = 1,
                 max_iter = 10, rho = 1.618, tau = 0.5,
                 e.rel = 1e-1, e.abs = 1e-1, phi = 1/200){
  n <- nrow(x.ori); p <- ncol(X)
  result <- admm_both(X, y, lambda = lamb, max_iter = max_iter, rho = rho, tau = tau,
                      e.rel = e.rel, e.abs = e.abs, phi = phi, M = M, G = G, K = K,
                      Args = NULL)
  residual <- abs(y - X %*% result$beta.hat)
  df <- sum(result$beta.hat != 0)
  BIC_cal <- log(sum(residual)) + df * log(log(n)) * log(p) / n
  list(hbic = BIC_cal, model = result)
}

#' Multiplier-bootstrap based lambda estimator (rank-style)
#' @param X design matrix
#' @param alpha0 1 - confidence level
#' @param c inflation factor
#' @param times number of draws
#' @export
est_lambda <- function(X, alpha0 = 0.1, c = 1.01, times = 1e2){
  dimn <- nrow(X)
  res <- NULL
  for (i in 1:times){
    epsilon_rank <- sample(1:dimn, dimn)
    xi <- 2*epsilon_rank - (dimn + 1)
    S <- (-2/dimn/(dimn-1)) * (t(X) %*% xi)
    res <- c(res, max(abs(S)))
  }
  as.numeric(stats::quantile(res, 1 - alpha0)) * c
}
