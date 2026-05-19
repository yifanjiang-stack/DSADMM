test_that("L1 and SCAD fits return finite coefficient vectors", {
  set.seed(11)
  beta <- c(1, 1, rep(0, 6))
  dat <- simulate_rank_example(n = 8, p = 8, rho = 0.3,
                               error = "normal", beta = beta)

  fit_l1 <- rank_admm_fit(
    dat$x, dat$y,
    penalty = "L1",
    lambda = 0.1,
    M = 2, G = 2, phi = 1 / 50,
    max_iter = 2,
    rank_transform = TRUE,
    standardize = TRUE
  )
  expect_length(fit_l1$beta, 8)
  expect_true(all(is.finite(fit_l1$beta)))

  fit_scad <- rank_admm_fit(
    dat$x, dat$y,
    penalty = "SCAD",
    lambda = 0.1,
    M = 2, G = 2, phi = 1 / 50,
    max_iter = 2,
    lla_max = 1,
    rank_transform = TRUE,
    standardize = TRUE
  )
  expect_length(fit_scad$beta, 8)
  expect_true(all(is.finite(fit_scad$beta)))
})

