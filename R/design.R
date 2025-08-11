#' Design matrix generators
#'
#' @param n number of rows
#' @param p number of columns
#' @param rho correlation parameter
#' @export
CS <- function(n, p, rho){
  X <- matrix(rnorm(n*p), n, p)
  X <- sqrt(1-rho)*X + sqrt(rho)*rnorm(n)
  X
}

#' @rdname CS
#' @export
AR <- function(n, p, rho){
  z <- matrix(NA_real_, n, p)
  for (j in 1:p) z[,j] <- rnorm(n)
  X <- matrix(NA_real_, n, p)
  X[,1] <- z[,1]
  for (j in 2:p){
    X[,j] <- rho*X[,j-1] + sqrt(1-rho^2)*z[,j]
  }
  X
}
