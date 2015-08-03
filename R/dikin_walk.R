#'Dikin Walk
#'
#'This function implements the Dikin Walk using the Hessian 
#'of the Log barrier function. Note that a $r$ of 1 guarantees
#'that the ellipsoid generated won't leave our polytope $K$ (see
#'Theorems online)
#'
#'@param A is the lhs of \eqn{Ax <= b}
#'@param b is the rhs of \eqn{Ax <= b}
#'@param x0 is the starting point
#'@param r is the radius of the ellipsoid (1 by default)
#'@param n is the number of points we want
#'@param analytic_center is a boolean indicating whether we want to 
#'       start at the analytic center of the polytope
#'
#'@return a number of sampled points that satisfy Ax <= b (matrix object, 
#'columns the points)
#'
#'@export


dikin_walk <- function(A, 
                       b, 
                       x0 = NULL, 
                       n, 
                       r = 1,
                       analytic_center = T) {
  
  ## assuming A is m * n matrix
  ## reference is here: http://mipt.ru/dcam/upload/06c/Seminar_MFTI-arpgzsi04lf.pdf
  ## Make the H(x) operator
  
  ## first, augment A | b 
  ## so that we could use "apply" for its speed
  A_b <- cbind (A, b)
  
  ## H(x) is the Hessian of the Log-barrier function at x
  
  H_x <- function(x) {
    
    ## making the D^2 matrix
    ## the lambda function is \frac{1}{b_i - a_i * x} , where a_i * x is the dot product of those two
    ## applying it to every row and then diag
    this.length <- ncol(A)
    D <- apply(A_b, 1, function(row) {1 / (tail(row, 1) - sum(row[1: (this.length)] * x))})
    
    ## diagonalize it
    ## still need to becareful of handling the case where 
    ## D = 0, so diag(0) returns an empty matrix 
    ## actually, probably can prove this: det(D) will never be zero in the case that 
    ## Ax <= b is bounded  ...?
    
    D_squared <- diag(D^2)
    
    ## this is the H(x) operator   
    
    return(t(A) %*% D_squared %*% A)  
    
  }
  
  ## D(x) is the diagonalized matrix of the log-barrier function of Ax <= b
  
  D_x <- function(x) {
    
    D <- apply(A_b, 1, function(row) {1 / (tail(row, 1) - sum(row[1: (length(row) - 1)] * x))})
    return(diag(D))
  }
  
  ## helper function:
  ## checks whether a point z is in Ellip(x)
  
  ellipsoid <- function(z, x) {
    
    return( sum(  (z-x) * ((H_x(x) %*% (z-x)))  ) <= r^2      )
    
  }
  
  ###### THE LINES ABOVE FINISH DEFINING THE ELLIPSOID
  ## now the sampling
  
  ## Given the H(x) operator, we must find points in the Dikin ellipsoid
  ## according to the equation
  ## Hd (lhs) = A^T * D * zeta (rhs), where zeta is a random point from the m-dimensional sphere
  ## length(b) should be m 
  
  result <- matrix(ncol = n+1, nrow = ncol(A))
  result[ , 1] <- x0
  current.point <- x0
  
  for (i in 1:n) {

    ## 1. Generate random point y in Ellip(x)
    ## KEY: MUST USE STANDARD NORMAL FUNCTION HERE
    ## If we use Unif[-1, 1], the sampling on the unit-sphere is not uniform
    ## and then the whole algorithm blows up 
    
    zeta <- rnorm(length(b), 0, 1)
    
    ## normalise to be on the m- unit sphere
    ## and then compute lhs as a m-vector
    
    zeta <- r * zeta / sqrt(sum(zeta * zeta))
    rhs <- t(A) %*% (D_x(current.point) %*% zeta )
    ## 
    
    #print(H_x(current.point))
    y <- solve(H_x(current.point), rhs)
    y <- y + current.point

   ## 2. Check whether x_0 is in Ellip(y)
   ## 3. Keep on trying y until condition satisfied
 
    while(!ellipsoid(current.point, y)) {
  #while(!dikin_ellipsoid(A = A, b = b, x0 = current.point, z = y, r = r)) {
     zeta <- rnorm(length(b), 0, 1)
      
      ## normalise to be on the m- unit sphere
      ## and then compute lhs as a m-vector
      zeta <- r * zeta / sqrt(sum(zeta * zeta))
      rhs <- t(A) %*% D_x(current.point) %*% zeta 
     
      ##
      #print(H_x(current.point))
      y <- solve(H_x(current.point), rhs)  ##THIS IS THE STEP THAT IS TAKING THE LONGEST
      y <- y + current.point

      if(ellipsoid(current.point, y)) {
      #if(dikin_ellipsoid(A = A, b = b, x0 = current.point, z = y, r =r )) {  

        ## det(A) / det(B) = 
        ## det(inv(B) %*% A)
        
        probability <- min(1, sqrt (det( solve(H_x(current.point)) %*% H_x(y))))
    
        bool <- sample(c(TRUE, FALSE), 1, prob = c(probability, 1-probability))
        if(bool) {
          
          ## get out of while loop
          ## break is dangerous, is there a better way of implementing this?
          
          break
        } 
      }
   }

  result[ , i+1] <- y
  current.point <- y
    

  }


  ## get rid of the columns which we don't sample
  ## should probably change this to check 
  
  cols <- which(!is.na(result[1,]))
  result <- result[,cols]
  return(result)

}
