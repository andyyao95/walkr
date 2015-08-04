#' Complete Solution
#' 
#' Given matrix equation \eqn{Aw=b}, find 
#' the basis representation of the 
#' complete solution (an affine transformation),
#' returning the homogeneous and particular solution
#' which describe the polytope in alpha-space
#' 
#' @param A is the lhs of the matrix equation \eqn{Ax=b}
#' @param b is the rhs of the matrix equation \eqn{Ax=b} 
#' @param randomize is the boolean indicating whether to 
#'        randomize the particular solution (i.e. randomize 
#'        obj. fxn with which it is generated)
#' 
#' @return a list object, with the first element($particular) as the particular solution
#'         and the the second element as a matrix with its columns containing
#'         the basis of the null space($homogeneous)
#'
#' @importFrom limSolve lsei
#' @importFrom MASS Null
#' 
#' @export
#'         

complete_solution <- function(A, b, randomize = FALSE) {
  
  
  ## Should do some checking here 
  
  ## must have an underdetermined matrix
  
  stopifnot(ncol(A) > nrow(A))
  
  
  ## Initialize the solution as a list 
  
  result <- list() 
  
  
  ## 1. Finding the particular solution 
  
  if(randomize == TRUE) {
    
    ## randomize objective function each time 
    ## (note: this actually doesnt guarantee a uniformly distribution of starting point, 
    ##  looking at the plots, they are very concentrated)
    
    ## these two lines are really just randomizing (to a certain extent): 
    ## the objective function to me minimized (which we don't care about)
    
    objfunc <- matrix(sample(c(-1, 1), ncol(A), replace = TRUE),
                      nrow = 1, ncol = ncol(A))
    const <- runif(1, -1, 1)
    
  }
  
  else {
    
    ## no randomization, just 1 1 1 ... 1 | 1 as the obj fxn
    
    objfunc <- matrix(rep(1, ncol(A)),
                      nrow = 1, ncol = ncol(A))
    const <- 1 
    
  } 
    
  ## particular solution is the 1st list object of the lsei fxn 
  ## here, we need to catch a warning that says the system has no solutions
  
  result$particular <- tryCatch (
    limSolve::lsei(A = objfunc, B = const, E = A, F = b)[[1]],
    warning = function(c) "not solvable"
  )
  
  ## solve warnings
  
  if(class(result$particular) == "character") {
    
    stop("The equality constraints cannot be all satisfied!
          Sampling in this sampling space is not possible!")
  }
  
  
  ## 2. Finding the homogeneous solution
  
  ## the Null function in MASS finds the null space of our solution
  
  result$homogeneous <- MASS::Null(t(A))
  
  ## Finally, the list object with particular and homogeneous solution is returned
  
  return(result)
  
  
}
