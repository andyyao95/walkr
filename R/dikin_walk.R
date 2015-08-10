#'Dikin Walk
#'
#'This function implements the Dikin Walk using the Hessian 
#'of the Log barrier function. Note that a $r$ of 1 guarantees
#'that the ellipsoid generated won't leave our polytope $K$ (see
#'Theorems online)
#'
#'@param A is the lhs of Ax <= b
#'@param b is the rhs of Ax <= b
#'@param x0 is the starting point (a list of points)
#'@param r is the radius of the ellipsoid (1 by default)
#'@param points is the number of points we want to sample
#'@param thin every thin-th point is stored
#'@param burn the first burn points are deleted
#'@param chains is the number of chains we run
#'
#' @return a list of chains of the sampled points, each chain
#'         being a matrix object with each column as a point
#'


dikin_walk <- function(A, 
                       b, 
                       x0 = list(), 
                       points, 
                       r = 1,
                       thin = 1,
                       burn = 0,
                       chains = 1) {
  
  stopifnot(points %% chains == 0)
  stopifnot(is.list(x0))
  #############################
  ## f stands for fast! 
  
  #1. rcppeigen_fprod(A, B) = A %*% b 
  #2. rcppeigen_fcrossprod(A, B) = t(A) %*% B
  #3. rcppeigen_ftcrossprod(A, B) = A %*% t(B)
  #4. rcppeigen_fsolve(A) = inverse of A   (solve(A))
  #5. rcppeigen_fdet(A) = determinant(A)
  
  ## first, augment A | b 
  A_b <- cbind (b, A)
  
  ## H(x) is the Hessian of the Log-barrier function of Ax <= b
  ## for more details, just google Dikin Walk
  
  H_x <- function(x) {
    
    ## making the D^2 matrix

    D <- as.vector(1/(A_b[,1] - rcppeigen_fprod(A_b[,-1], x)))

    ## t(A) %*% (D^2 %*% A)
    
    return(rcppeigen_fcrossprod(A, rcppeigen_fprod(diag(D^2), A)))
    
  } 
  
  ## D(x) is the diagonalized matrix of 1 over the log-barrier function of Ax <= b
  
  D_x <- function(x) {

    return(diag(as.vector(1/(A_b[,1] - rcppeigen_fprod(A_b[,-1], x)))))
  } 
  
  ## checks whether a point z is in Dikin Ellip centered at x
  
  ellipsoid <- function(z, x) {
    
    ## as.numeric converts the expression into an atom, so we get boolean
    ## it's just checking (z-x)^T %*% H_x %*% (z-x) <= r^2  
  
    return( as.numeric(rcppeigen_fcrossprod(z-x, rcppeigen_fprod(H_x(x), (z-x)))) <= r^2)
    
  } 
  
  ###### THE LINES ABOVE FINISH DEFINING THE ELLIPSOID
  ## now the sampling
  
  ## initialize return matrix
  ## set the starting point as the current point
  answer <- list()
  
  ##total points for each indiv chain
  
  for (j in 1:chains) {
  
    total.points <- ( (points / chains)  * thin + burn) 
  
    
    result <- matrix(ncol = total.points, nrow = ncol(A))
    result[ , 1] <- x0[[j]]
    current.point <- x0[[j]]
    this.length <- length(b)
    
    for (i in 2:total.points) {
      
      ## 1. Generate random point y in Ellip(x)
      ## KEY: MUST USE STANDARD NORMAL FUNCTION HERE, RUNIF IS NOT UNIFORM AFTER TRANSFORMATION
      ## see vignette for details
      
      zeta <- stats::rnorm(this.length, 0, 1)
      
      ## normalise to be on the m- unit sphere
      ## and then compute lhs as a m-vector
      
      ## essentially: Hd = t(A) %*% D^2 %*% zeta
      ## solving for d gives us a uniformly random vector in the ellipsoid centered at x 
      ## the y = x_0 + d is the new point 
      
      zeta <- r * zeta / sqrt(as.numeric(rcppeigen_fcrossprod(zeta,zeta)))
      rhs <- rcppeigen_fcrossprod(A, rcppeigen_fprod(D_x(current.point), zeta))
      
      y <- rcppeigen_fprod(rcppeigen_fsolve(H_x(current.point)), rhs) + current.point 
      
  
      ## 2. Check whether x_0 is in Ellip(y)
      ## 3. Keep on trying y until condition satisfied
      
      while(!ellipsoid(current.point, y)) {
        
        ## exact same set of procedures as above
        
        zeta <- stats::rnorm(this.length, 0, 1)
        zeta <- r * zeta / sqrt(sum(zeta * zeta))
        rhs <- rcppeigen_fcrossprod(A, rcppeigen_fprod(D_x(current.point), zeta))
        y <- rcppeigen_fprod(rcppeigen_fsolve(H_x(current.point)), rhs) + current.point 
        
        if(ellipsoid(current.point, y)) {
          
          ## det(A)/det(B) = det(B^-1 A)
          ## acceptance rate according to probability formula. see paper for detail
          
          probability <- min(1, sqrt (rcppeigen_fdet( rcppeigen_fprod(
            rcppeigen_fsolve(H_x(current.point)),H_x(y)))))
          
          bool <- sample(c(TRUE, FALSE), 1, prob = c(probability, 1-probability))
          
          if(bool) {
            
            ## perhaps, there is a better way to handle break here?
            
            break
          } 
          
        }
      }
      
      ## appending on the result
      
      result[ , i] <- y
      current.point <- y
      
      
    }
    
    ## NEED TO HANDLE THE CASE WHEN ALPHA IS JUST 1 DIMENSIONAL
  
    if(dim(result)[1] == 1) {
      
      ## first, delete out the number of points that we want to burn
      ## second, only take every thin-th point
      
      result <- matrix(result[, (1+burn) : total.points], nrow = 1)
      result <- matrix(result[ , (1:(points/chains))*thin], nrow = 1)
    }
    
    else {
      
      
      ## first, delete out the number of points that we want to burn
      ## second, only take every thin-th point
      
      result <- result[, (1+burn) : total.points]
    
      result <- result[ , (1:(points/chains))*thin]
    }
    
    ## appending on 1 chain onto the result
    
    answer[[j]] <- result
  }
  
  return(answer)
}
