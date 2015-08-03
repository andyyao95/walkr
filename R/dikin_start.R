#' Dikin Start
#' 
#' Given \eqn{Ax <= b}, which defines a convex
#' polytope, this function picks n random
#' starting "center" points using linear programming. 
#' 
#' @param A is the lhs of \eqn{Ax <= b}
#' @param b is the rhs of \eqn{Ax <= b}
#' @param n is the number of points we want to return
#' @param average is the number of boundary points we want 
#'        to take the average of
#' 
#' @return a matrix, with each column as a point
#' 
#' @export


dikin_start <- function(A, 
                        b, 
                        n = 1, 
                        average = 10) {
  
  ## initialize the return matrix, each column is 1 point
  ## the points have dimension equal to the rows of A
  
  result <- matrix(ncol = n, nrow = ncol(A))
 
  ## creating every point
  
  for (i in 1:n) {
  
    ## initialize local variable to store boundary points
    ## so that we could take an average in the end
    
    new_x0 <- numeric()
  
    ## average is the number of boundary points we want to take the average of 
    ## the higher this number is, the more likely our n points will be closer to each other
    ## 
    
    for(j in 1:average) {
        
        ## these two lines randomize 
        
        objfunc <- matrix(sample(runif(1,-1, 1), ncol(A), replace = TRUE),
                          nrow = 1, ncol = ncol(A))
        const <- runif(1, -1, 1)
        
        ## suppresswarnings because we don't specific a equality constraint,
        ## we only care about Ax <= b 
        
        ## in terms of the lsei, it actually takes in Gx >= h
        ## thus, we pass in -Ax >= -b , which is equivalent to Ax <= b
        
        new_x0 <- rbind(new_x0, tryCatch (    
          suppressWarnings(
            limSolve::lsei(A = objfunc, B = const, G = -1* A, H = -1 * b + 0.000001)[[1]]),
          error = function(c) stop("The inequality constraints cannot be all satisfied!
                                   Sampling in this solution space is not possible!")
      ))
    }
    
    ## take the average of those boundary points ==> to obtain the center
    
    result[, i] <- apply(new_x0,  2, mean)
 
  }
  
  
  return(result)
}