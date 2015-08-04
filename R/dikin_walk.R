#'Dikin Walk
#'
#'This function implements the Dikin Walk when given a
#'convex polytope defined by \eqn{Ax \le b} and a starting 
#'point x0.
#'Note that a $r$ of 1 guarantees
#'that the ellipsoid generated won't leave our polytope $K$ (see
#'Theorems online)
#'
#'@param A is the lhs of \eqn{Ax \le b}
#'@param b is the rhs of \eqn{Ax \le b}
#'@param x0 is the starting point
#'@param r is the radius of the Dikin ellipsoid (1 by default)
#'@param n is the number of points we want to sample
#'
#'@return a number of sampled points that satisfy Ax <= b (matrix object, 
#'columns the points)
#'


dikin_walk <- function(A, 
                       b, 
                       x0,  
                       n, 
                       r = 1) {
  
  ## assuming A is m * n matrix
  ## reference is here: http://mipt.ru/dcam/upload/06c/Seminar_MFTI-arpgzsi04lf.pdf

  ## initialize return matrix
  ## set the starting point as the current point
  
  result <- matrix(ncol = n+1, nrow = ncol(A))
  result[ , 1] <- x0
  current.point <- x0
  
  for (i in 1:n) {

    ## 1. Generate random point y in Ellip(x)
    ## KEY: MUST USE STANDARD NORMAL FUNCTION HERE, RUNIF IS NOT UNIFORM AFTER TRANSFORMATION
    ## see vignette for details
    
    zeta <- rnorm(length(b), 0, 1)
    
    ## normalise to be on the m- unit sphere
    ## and then compute lhs as a m-vector
    
    ## essentially: Hd = t(A) %*% D^2 %*% zeta
    ## solving for d gives us a uniformly random vector in the ellipsoid centered at x 
    ## the y = x_0 + d is the new point 
    
    zeta <- r * zeta / sqrt(sum(zeta * zeta))
    
    rhs <- t(A) %*% (diag_logbarrier(A = A, b = b, x = current.point) %*% zeta )
     
    y <- solve(hessian_logbarrier(A = A, b = b, x = current.point), rhs)
    y <- y + current.point

    ## 2. Check whether x_0 is in Ellip(y)
    ## 3. Keep on trying y until condition satisfied
 

    while(!dikin_ellipsoid(A = A, b = b, x0 = y, z = current.point, r = r)) {
      
      zeta <- rnorm(length(b), 0, 1)
      
      ## exact same set of procedures as above
     
      zeta <- r * zeta / sqrt(sum(zeta * zeta))
      rhs <- t(A) %*% (diag_logbarrier(A = A, b = b, x = current.point) %*% zeta )
      y <- solve(hessian_logbarrier(A = A, b = b, x = current.point), rhs)
      y <- y + current.point

      if(dikin_ellipsoid(A = A, b = b, x0 = y, z = current.point, r = r)) {  

        ## det(A)/det(B) = det(B^-1 A)
        ## acceptance rate according to probability formula. see paper for detail
        
        probability <- min(1, sqrt (det( solve(hessian_logbarrier(A = A,
                                                                  b = b, x = current.point)) %*% 
                                                                  hessian_logbarrier(A = A, b = b, x = y))))
    
        bool <- sample(c(TRUE, FALSE), 1, prob = c(probability, 1-probability))
        
        if(bool) {
          
          ## perhaps, there is a better way to handle break here?
          
          break
        } 
      }
   }
   
   ## append result
  
   result[ , i+1] <- y
   current.point <- y
    
  }


  ## get rid of the columns which we don't sample (or fail to sample?)
  ## for safety
  
  cols <- which(!is.na(result[1,]))
  result <- result[,cols]
  
  return(result)

}
