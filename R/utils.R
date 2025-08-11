#' Proximal and utility helpers
#'
#' @param u,a numeric vectors
#' @param xi,alpha numeric vectors
#' @param quan quantile level in (0,1)
#' @export
prox_beta <- function(u,a){
  sign(u)*pmax(0, (abs(u)-a))
}

#' @rdname prox_beta
#' @export
prox_z <- function(xi,alpha,quan){
  pmax(xi-quan*alpha,0) - pmax(-xi-(1-quan)*alpha,0)
}

#' Euclidean norm
#' @param x numeric vector
#' @export
l2norm <- function(x){
  sqrt(sum(x^2))
}

#' L1 norm
#' @param x numeric vector
#' @export
l1norm <- function(x){
  sum(abs(x))
}

#' Matrix splitter (internal utility)
#' @param M matrix to split
#' @param r number of row blocks
#' @param c number of column blocks
#' @return A nested list indexed as [row_block, col_block]
#' @export
matsplitter <- function(M, r, c) {
  stopifnot(is.matrix(M), r >= 1, c >= 1)
  splitMatrix <- function(mat, nrow) {
    split.data.frame(t(mat), ceiling(1:ncol(mat)/ncol(mat)*nrow))
  }
  sapply(splitMatrix(M, c), splitMatrix, r)
}
