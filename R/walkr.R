#' walkr 
#' 
#' Given Ax = b, sample points from the intersection of Ax = b 
#' with the n-simplex (\eqn{\sum(x) == 1}, \eqn{x_i >= 0}). The
#' current MCMC sampling methods supported are "hit-and-run" and 
#' "dikin"
#' 
#' @param A is the lhs of the matrix equation A
#' @param b is the rhs of the matrix equation b
#' @param n is the number of points we want to sample
#' @param method is the MCMC sampling method. Please enter "hit-and-run" or "dikin"
#' 
#' @return a list object. The first element is a matrix with its columns as the sampled
#'         points. The second element is the diagnostics of the MCMC sampling. 
#'         
#' @export 
#' 

walkr <- function(A, b, n, method) {
  
  ## 0. Should do some checking here
  
  ## 1. regardless of method, we need to perform the affine transformation which
  ## takes us from x-space (Ax = b) into the alpha-space in which the polytope
  ## described is Ax <= b. From there, we could perform the sampling
  
  z <- complete_solution(A,b, randomize = T)
  
  ## break it up into null space(homogeneous) and particular solution
  
  particular  <- z$particular
  homogeneous <- z$homogeneous
  
  ## Homogeneous %*% alpha >= -vp
  ## -Homogeneous %*% alpha <= vp    (Ax <= b form)
  
  new_A <- -homogeneous
  new_b <- particular
  
  
  ## 2. Find starting point within convex polytope
  
  ## Finding a starting point to run our chain
  
  x0 <- start_point(A = new_A, b = new_b, n = 1, average = 20)
  start.point <- as.vector(x0)
  
  
  ## 3. The sampling
  
  
  if(method == "dikin") {
    
    alphas <- dikin_walk(A = new_A, b = new_b, n = n-1, r = 1, x0 = start.point)
    
    answer <- apply(alphas, 2, function(x) { homogeneous %*% x + particular  })
    
    return(answer)
  }
  
  else if (method == "hit-and-run") {
    
    stop("not yet implemented!")
  }
  
  
  else{
    stop("Sampling method must be \"hitandrun\" or \"dikin\".")
  }
  
  return(0)
  
}