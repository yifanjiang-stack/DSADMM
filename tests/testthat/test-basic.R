test_that("creating_x/y sizes are correct", {
  x <- matrix(1:9, nrow = 3)
  y <- 1:3
  expect_equal(nrow(creating_x(x)), 3)   # 3*2/2
  expect_equal(ncol(creating_x(x)), ncol(x))
  expect_equal(length(creating_y(y)), 3) # 3*2/2
})

test_that("AR generator returns correct dims", {
  X <- AR(10, 5, rho = 0.3)
  expect_true(is.matrix(X))
  expect_equal(dim(X), c(10,5))
})
