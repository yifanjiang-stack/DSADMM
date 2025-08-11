# DSADMM

Helpers and an ADMM solver for quantile regression (double-splitting style).

## Install

```r
# install.packages("remotes")
remotes::install_github("YOUR_GITHUB/DSADMM")
```

## Usage

```r
library(DSADMM)

# toy design
set.seed(1)
n <- 200; p <- 50
X <- AR(n, p, rho = 0.5)
beta_true <- numeric(p); beta_true[c(6, 10, 25)] <- c(1, 0.5, -0.8)
y <- as.vector(X %*% beta_true + rnorm(n))

fit <- admm_both(X, y, lambda = 0.1, max_iter = 100, rho = 1.618, tau = 0.5,
                 e.rel = 1e-3, e.abs = 1e-3, phi = 1/200, M = 2, G = 1, K = 1,
                 beta_true = beta_true, idx_true = which(beta_true != 0))

str(fit)
```
