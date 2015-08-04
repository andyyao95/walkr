#' diagonalized log-barrier 
#' 
#' Computes the diagonal matrix from the vector obtained 
#' from 1 / the log-barrier function
#' 
#'@param A is the lhs of \eqn{Ax <= b}
#'@param b is the rhs of \eqn{Ax <= b}
#'@param x is the starting point
#'
#'@return a diagonalized matrix of 1/logbarrier function
#'

diag_logbarrier <- function(A, b, x) {
  
  ## augment A and b so we can use apply
  
  A_b <- cbind(A, b)
  
  ## each element in the vector is 1 / (b_i - a_i^T x)
  
  D <- apply(A_b, 1, function(row) {1 / (tail(row, 1) - sum(row[1: (length(row) - 1)] * x))})
  
  return(diag(D))
  
}