#'Dikin Walk
#'
#'This function implements the Dikin Walk using the Hessian 
#'of the Log barrier function. Note that a $r$ of 1 guarantees
#'that the ellipsoid generated won't leave our polytope $K$ (see
#'Theorems online)
#'
#'@param A is the lhs of Ax <= b
#'@param b is the rhs of Ax <= b
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


Dikin_Walk <- function(A, 
                       b, 
                       x0 = NULL, 
                       n, 
                       r = 1,
                       analytic_center = T) {
  
  ## first, augment A | b 
  A_b <- cbind (b, A)
  
  ## H(x) is the Hessian of the Log-barrier function of Ax <= b
  ## for more details, just google Dikin Walk
  
  H_x <- function(x) {
    
    ## making the D^2 matrix
    ## the lambda function is \frac{1}{b_i - a_i * x} , where a_i * x is the dot product of those two
    ## applying it to every row and then diag

    D <- as.vector(1/(A_b[,1] - fprod(A_b[,-1], x)))

    ## t(A) %*% (D^2 %*% A)
    
    return(fcrossprod(A, fprod(diag(D^2), A)))
    
  } ##DONE
  
  ## D(x) is the diagonalized matrix of the log-barrier function of Ax <= b
  
  D_x <- function(x) {

    #D <- as.vector(1/(A_b[,1] - fprod(A_b[,-1], x)))
    return(diag(as.vector(1/(A_b[,1] - fprod(A_b[,-1], x)))))
  } 
  
  ## helper function:
  ## checks whether a point z is in Ellip(x)
  
  ellipsoid <- function(z, x) {
    
    ## as.numeric converts the expression into an atom, so we get boolean
    return( as.numeric(fcrossprod(z-x, fprod(H_x(x), (z-x)))) <= r^2)
    
  } 
  
  ###### THE LINES ABOVE FINISH DEFINING THE ELLIPSOID
  ## now the sampling
  result <- matrix(ncol = n+1, nrow = ncol(A))
  result[ , 1] <- x0
  current.point <- x0
  this.length <- length(b)
  
  for (i in 1:n) {
    if(i %% 100 == 0) {print(i)}
    ## 1. Generate random point y in Ellip(x)
    ## KEY: MUST USE STANDARD NORMAL FUNCTION HERE, I DONT KNOW WHY BUT RUNIF BLOWS THINGS UP
    
    zeta <- rnorm(this.length, 0, 1)
    
    ## normalise to be on the m- unit sphere
    ## and then compute lhs as a m-vector
    
    zeta <- r * zeta / sqrt(as.numeric(fcrossprod(zeta,zeta)))
    rhs <- fcrossprod(A, fprod(D_x(current.point), zeta))
    #print(H_x(current.point))
    y <- fprod(fsolve(H_x(current.point)), rhs) + current.point 
    

    ## 2. Check whether x_0 is in Ellip(y)
    ## 3. Keep on trying y until condition satisfied
    
    while(!ellipsoid(current.point, y)) {
      #print("inwhile")
      zeta <- rnorm(this.length, 0, 1)
      zeta <- r * zeta / sqrt(sum(zeta * zeta))
      rhs <- fcrossprod(A, fprod(D_x(current.point), zeta))
      #print(H_x(current.point))
      y <- fprod(fsolve(H_x(current.point)), rhs) + current.point 
      
      if(ellipsoid(current.point, y)) {
        
        ## genius.... det(A)/det(B) = det(B^-1 A)
        
        probability <- min(1, sqrt (fdet( fprod(fsolve(H_x(current.point)),H_x(y)))))
        
        bool <- sample(c(TRUE, FALSE), 1, prob = c(probability, 1-probability))
        
        if(bool) {
          break
        } 
        
      }
    }
    
    result[ , i+1] <- y
    current.point <- y
    
    
  }
  
  
  ## get rid of the columns which we don't sample
  
  #cols <- which(!is.na(result[1,]))
  #result <- result[,cols]
  return(result)
  
}
