#' Pairwise difference constructors
#'
#' `creating_x()` builds all pairwise row differences of a matrix.
#' `creating_y()` builds all pairwise differences of a vector.
#'
#' @param x A numeric matrix with \eqn{n} rows and \eqn{p} columns.
#' @param y A numeric vector of length \eqn{n}.
#' @return
#' `creating_x()` returns a matrix with \eqn{n(n-1)/2} rows and \eqn{p} columns.
#' `creating_y()` returns a vector of length \eqn{n(n-1)/2}.
#' @examples
#' x <- matrix(1:6, nrow = 3)
#' creating_x(x)
#' y <- c(1, 4, 9)
#' creating_y(y)
#' @export
creating_x <- function(x){
  stopifnot(is.matrix(x))
  x_tilde <- matrix(ncol = ncol(x), nrow = nrow(x) * (nrow(x) - 1) / 2)
  k <- 1L
  for (i in 1:(nrow(x)-1)){
    for (j in (i+1):nrow(x)){
      x_tilde[k,] <- x[i,]-x[j,]
      k <- k + 1L
    }
  }
  x_tilde
}

#' @rdname creating_x
#' @export
creating_y <- function(y){
  stopifnot(is.numeric(y))
  y <- as.vector(y)
  y_tilde <- vector(mode = storage.mode(y), length = length(y)*(length(y)-1)/2)
  k <- 1L
  for (i in 1:(length(y)-1)){
    for (j in (i+1):length(y)){
      y_tilde[k] <- y[i]-y[j]
      k <- k + 1L
    }
  }
  y_tilde
}

# Optional byte-code compilation
#' @importFrom compiler cmpfun
creating_x <- compiler::cmpfun(creating_x)
#' @importFrom compiler cmpfun
creating_y <- compiler::cmpfun(creating_y)
