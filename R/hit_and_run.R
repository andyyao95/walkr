#' Hit and Run 
#' 
#' This function provides a wrapper 
#' for the har function of the 
#' hit-and-run package
#' 
#' @param A is the lhs of the matrix equation A
#' @param b is the rhs of the matrix equation b
#' @param points is the number of points we want to sample
#' @param thin -- every thin-th point is taken
#' @param burn the first burn number of points are deleted
#' 
#'@return a number of sampled points that satisfy Ax <= b (matrix object, 
#'        columns the points)
#' 
#' @importFrom hitandrun har 


hit_and_run <- function(A, 
                        b, 
                        x0, 
                        points, 
                        thin = 1,
                        burn = 0) {
  
  
  ## need some checking here
  
  ## initialize the list of constraints
  ## should also do some checking
  
  constr <- list(constr = A, rhs = b, dir = rep("<=", nrow(A)))
  
  ## getting samples
  
  total.points <- ( points + burn ) * thin 
  
  result <- t(hitandrun::har(x0, constr, N = total.points, 
                   thin = thin, )$samples)
  
  ## NEED TO HANDLE THE CASE WHEN ALPHA IS JUST 1 DIMENSIONAL
  
  if(dim(result)[1] == 1) {
    
    result <- matrix(result[, (1+burn) : ncol(result)], nrow = 1)
    
  }
  
  else {
    result <- result[ , (1+burn) : ncol(result)]
    
  }
  
  return(result)
  
  
} 