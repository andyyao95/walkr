#' walkr 
#' 
#' Given Ax = b, sample points from the intersection of Ax = b 
#' with the n-simplex (\eqn{\sum x = 1}, \eqn{x_i \ge 0}). The
#' current MCMC sampling methods supported are "hit-and-run" and 
#' "dikin"
#' 
#' @param A is the lhs of the matrix equation A
#' @param b is the rhs of the matrix equation b
#' @param n is the number of points we want to sample
#' @param method is the MCMC sampling method. Please enter "hit-and-run" or "dikin"
#' 
#' @return A matrix with its columns as the sampled
#'         points. 
#'         
#' @importFrom hitandrun har          
#' @export 
#' 

walkr <- function(A, 
                  b, 
                  n, 
                  method = "dikin") {
  
  ## 0. Should do some checking here
  
  ## 1. regardless of method, we need to perform the affine transformation which
  ## takes us from x-space (Ax = b) into the alpha-space in which the polytope
  ## described is Ax <= b (denoted below as new_A, new_b)
  ## From there, we could perform the sampling
  
  ## break it up into null space(homogeneous) and particular solution
  
  z <- complete_solution(A,b, randomize = T)
  particular  <- z$particular
  homogeneous <- z$homogeneous
  
  ## Homogeneous %*% alpha >= -vp
  ## -Homogeneous %*% alpha <= vp    (Ax <= b form)
  
  new_A <- -homogeneous
  new_b <- particular
  
  
  ## 2. Find starting point within convex polytope
  
  x0 <- start_point(A = new_A, b = new_b, n = 1, average = 20)
  start.point <- as.vector(x0)
  
  
  ## 3. The sampling
  
  
  if(method == "dikin") {
    
    ## sampling in alpha space
    ## n = n - 1 because dikin takes starting point as the 1st sampled point
    
    alphas <- dikin_walk(A = new_A, b = new_b, n = n-1, r = 1, x0 = start.point)
    
    ## convert back into x-space
    
    answer <- apply(alphas, 2, function(x) { homogeneous %*% x + particular  })
    
    return(answer)
  }
  
  else if (method == "hit-and-run") {
    
    ## make the constraints in the format that the hitandrun package wants
    
    constr <- list(constr = new_A, rhs = new_b, dir = rep("<=", nrow(new_A)))
     
    ##again, sampling alphas
    
    alphas <- t(hitandrun::har(start.point, constr, N = n, 
                              thin = 1, )$samples)
    answer <- apply(alphas, 2, function(x) { homogeneous %*% x + particular  })
    
    return(answer)
  }
  
  
  else{
    stop("Sampling method must be \"hitandrun\" or \"dikin\".")
  }
  
  ## for safety. function should never hit here
  
  return(0)
  
}