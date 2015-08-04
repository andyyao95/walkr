#'Hessian of Log Barrier
#'
#'At point x0, this function computes the Hessian
#'of the log barrier function for the 
#'convex polytope that is defined by \eqn{Ax \le b}.
#'
#'@param A is the lhs of \eqn{Ax \le b}
#'@param b is the rhs of \eqn{Ax \le b}
#'@param x is the current point we're at
#'
#'@return the Hessian matrix for the specific parameters

hessian_logbarrier <- function(A, b, x) {
  
  ## first, augment A | b 
  ## so that we could use "apply" for its speed
  A_b <- cbind (A, b)
  
  ## making the D^2 matrix
  ## the lambda function is \frac{1}{b_i - a_i * x} , where a_i * x is the dot product of those two
  ## applying it to every row and then diag
  this.length <- ncol(A)
  
  D <- apply(A_b, 1, function(row) {1 / (tail(row, 1) - sum(row[1: (this.length)] * x))})
  
  ## diagonalize it
  ## still need to becareful of handling the case where 
  ## D = 0, so diag(0) returns an empty matrix 
  
  D_squared <- diag(D^2)
  
  ## this is the H(x) operator   
  
  return(crossprod(A, (D_squared %*% A) ))
    
}