Dikin_single_chain <- function(x0,
                               total.points,
                               A,
                               b,
                               r,
                               A_b,
                               burn,
                               chains,
                               points,
                               thin,
                               c){
  
  
  H_extra <- function (c, x) {
    # function returns Hessian of exp(X-C)
    x <- as.vector(x)
    distance <- sum((x - c)^2)
    return(4*exp(distance)*(x - c) %*% t(x - c) +
             2*exp(distance)*diag(length(x)))  
  }
  
  H_x <- function(x, c) {
    
    ## making the D^2 matrix
    
    D <- as.vector(1/(A_b[,1] - rcppeigen_fprod(A_b[,-1], x)))
    
    ## t(A) %*% (D^2 %*% A)
      
    return(rcppeigen_fcrossprod(A, rcppeigen_fprod(diag(D^2), A)) + 
             H_extra(c, x))
    
  } 
  
  ## D(x) is the diagonalized matrix of 1 over the log-barrier function of Ax <= b
  
  D_x <- function(x) {
    
    return(diag(as.vector(1/(A_b[,1] - rcppeigen_fprod(A_b[,-1], x)))))
  } 
  
  ## checks whether a point z is in Dikin Ellip centered at x
  
  ellipsoid <- function(z, x) {
    
    ## as.numeric converts the expression into an atom, so we get boolean
    ## it's just checking (z-x)^T %*% H_x %*% (z-x) <= r^2  
    
    #     \subsection{How to pick a random point uniformly from a Dikin Ellipsoid?}
    #     
    #     Let's say, we now have $D_{x}^r$, the Dikin Ellipsoid centered at $x$ with radius $r$. 
    #     
    #     \begin{enumerate}
    #     
    #     \item{generate $\zeta$ from the $n$ dimensional Standard Gaussian (i.e. \texttt{zeta = rnorm(n,0,1)})}
    #     \item{normalize $\zeta$ to be on the $n$ dimensional ball with radius $r$, that is:}
    #     \subitem{$\zeta \quad = \quad <x_1,x_2,...,x_n> \quad \rightarrow \quad <\frac{rx_1}
    #     {\sqrt{x_1^2+x_2^2+...+x_n^2}}, 
    #     \frac{rx_2}{\sqrt{x_1^2+x_2^2+...+x_n^2}}, ...... , 
    #     \frac{rx_n}{\sqrt{x_1^2+x_2^2+...+x_n^2}}>$}
    #     \item{Solve for $d$ in the matrix equation $H_x d = A^TD\zeta$ (note, as long as $x_0$ is not on the boundary
    #     of our polytope $K$, $H_x$ will be non-singular, thus, $d$ will always be unique)}
    #     \item{$y = x_0 + d$ is our randomly sampled point from $D_x^r$}
    #     
    #     \end{enumerate}
    #     
    #     
    #     1) will clean this up later
    #     2) should move these 3 functions into separate files with documentation...
    #     
    
    return( as.numeric(rcppeigen_fcrossprod(z-x, rcppeigen_fprod(H_x(x, c), (z-x)))) <= r^2)
    
  } 
  
  
  result <- matrix(ncol = total.points, nrow = ncol(A))
  result[ , 1] <- x0
  current.point <- x0
  this.length <- length(b)
  
  for (i in 2:total.points) {
    
    ## 2. Check whether x_0 is in Ellip(y)
    ## 3. Keep on trying y until condition satisfied
    
    bool <- FALSE
    
    ## will always go into the while-loop for the first time, due to lazy evaluation
    
    while(!bool || !ellipsoid(current.point, y)) {
      
      ## exact same set of procedures as above
      
      zeta <- stats::rnorm(this.length, 0, 1)
      zeta <- r * zeta / sqrt(sum(zeta * zeta))
      rhs <- rcppeigen_fcrossprod(A, rcppeigen_fprod(D_x(current.point), zeta))
      
      # declare a new variable to save one computation later
      inverseTemp <- rcppeigen_fsolve(H_x(current.point, c))
      
      y <- rcppeigen_fprod(inverseTemp, rhs) + current.point 
      
      if(ellipsoid(current.point, y)) {
        
        ## det(A)/det(B) = det(B^-1 A)
        ## acceptance rate according to probability formula. see paper for detail
        
        probability <- min(1, sqrt(rcppeigen_fdet(rcppeigen_fprod(inverseTemp, H_x(y, c)))))
        
        bool <- sample(c(TRUE, FALSE), 1, prob = c(probability, 1 - probability))
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
    
    ## we take the floor function because we took the ceiling above
    ## so in the case that multiplying by burn doesn't result in an integer
    ## we returning the correct number of points
    
    result <- matrix(result[, (floor(burn*total.points)+1) : total.points], nrow = 1)
    result <- matrix(result[ , (1:(points/chains))*thin], nrow = 1)
  }
  
  else {
    
    
    ## first, delete out the number of points that we want to burn
    ## second, only take every thin-th point
    
    ## same as above
    
    result <- result[, (floor(burn*total.points)+1) : total.points]
    
    result <- result[ , (1:(points/chains))*thin]
  }
  
  return(result)

}
