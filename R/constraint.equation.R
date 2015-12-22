#' Creates a constraint equation based on the input data frame.
#' 
#' @param x data frame containing needed input data
#' @param weight.var character name of the column of the input weights
#' @param match.var character vector of names of columns of 'data' we wish to 
#'   match on
#' @param replace logical indicating whether or not bservations weighted in the 
#'   original weight.var are allowed positive weight in the output.
#'   
#' @return A list with two named components: A and b, representing the
#'   components of the constraint equation \eqn{Ax = b}
#'   
#' @author David Kane \email{<dave.kane@@gmail.com>}
#' @export
#' 
constraint.equation <- function(x, weight.var, match.var, replace){
  
  dummy <- function(vec){
    
    ## sort in order to match up match up exposures when printing
    names <- sort(unique(vec))
    mat <- matrix(rep(0, length(vec)*length(names)), nrow = length(names))
    for(i in 1:nrow(mat)) {
      mat[i,][which(vec == names[i])] <- 1
    }
    return(mat)
  }
  
  
  
  ## Intialize list that will be turned into a matrix with do.call
  
  Alist = list()
  
  ## Include the continuous and discrete variables in Alist. Need to be careful
  ## in deciding just what a 0/1 variable means, for example.
  
  for(i in 1:length(match.var)){
    if(is.numeric(x[[match.var[i]]])){
      ## continuous
      Alist[[i]] = x[[match.var[i]]]
      
    } else {
      ## discrete
      Alist[[i]] = dummy(x[[match.var[i]]])
    }
  }
  ## do.call on Alist makes constraint matrix A
  
  A = do.call(rbind, Alist)
  
  ## b is the constraint matrix
  b = A %*% x[[weight.var]]
  
  ## attach the "match sum" constraint, redundanies no longer matter
  sumlimit = sum(x[[weight.var]])
  A = rbind(A, rep(1, ncol(A)))
  b = c(b, sumlimit)
  
  if(!replace) {
    if(sum(x[[weight.var]] > 0) >= nrow(x) ) stop("All rows are weighted, set replace = TRUE.")
    ## remove columns corresponding to variables that have weight in the original 
    A = A[,-which(x[[weight.var]] > 0)]
  }
  return(list(A = A, b = b)) 
}