#' Dikin Ellipsoid
#' 
#' Given specific parameters,
#' this function checks whether a point
#' z is in the specified Dikin Ellipsoid 
#' centered at x0
#' 
#'@param A is the lhs of \eqn{Ax <= b}
#'@param b is the rhs of \eqn{Ax <= b}
#'@param x0 is the starting point
#'@param r is the radius of the ellipsoid 
#'@param z is the point we want to test
#'
#'@return a boolean indicating whether z is in the ellipsoid
#'
#'@export 
#'

dikin_ellipsoid <- function(A, b, x0, r, z) {
    
    ## The Dikin Ellipsoid centered at x0 with radius r is defined as:
    ## All the points z which satisfy
  
    ## (z-x)^T H_{x0} (z-x) <= r^2
    ## where H_{x0} is the Hessian of the log-barrier function at x_0 for Ax <= b 
  
    return( sum(  (z-x0) * ((hessian_logbarrier(A = A, b = b, x = x0) %*% (z-x0)))  ) <= r^2      )
  
}