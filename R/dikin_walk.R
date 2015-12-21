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
#' @examples
#' \dontrun{
#' ## note that this Ax <= b is different from Ax=b that the 
#' ## user specifies for walkr (see transformation section in vignette) 
#' dikin_walk(A = A, b = b, x0, points = 100, 
#'            r = 1thin = 1, burn = 0, chains = 1) 
#' } 

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
    
    return( as.numeric(rcppeigen_fcrossprod(z-x, rcppeigen_fprod(H_x(x), (z-x)))) <= r^2)
    
  } 
  
  ###### THE LINES ABOVE FINISH DEFINING THE ELLIPSOID
  ## now the sampling
  
  ## initialize return matrix
  ## set the starting point as the current point
  answer <- list()
  
  ##total points for each indiv chain
  
  total.points <- ceiling( (points / chains)  * thin * (1/(1-burn))) 
  
  for (j in 1:chains) {
  
    ## total points is : points * thin * 1/(1-burn) / chains 
    ## because burn-in is a percentage, we must take the CEILING function
    ## to sample more than we need (in the case where dividing by 1-burn
    ## does not return an integer)
    answer[[j]] <- Dikin_single_chain(total.points, A, b, x0[[j]],
                                      r, A_b, burn, chains, points, thin)
    
    ## initializing the return matrix 
    
    
    ## appending on 1 chain onto the result
    
  }
  
  return(answer)
}
