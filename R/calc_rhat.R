#' R hat 
#' 
#' This function calculates the rhat
#' of each parameter given the list of chains 
#' that walkr produces. Since this is just an 
#' internal function, I'll document it more later.
#' 
#' @param x is the list of chains 
#' 
#' @return a vector of rhats 


calc_rhat <- function(x) {
  
  ## first, since the test cases for walkr passed
  ## we can safely assume here that x is of the correct format
  ## which we require it to be
  
  ## m is the number of chains
  
  m <- length(x)
  
  ## n is the number of points in each chain
  
  n <- dim(x$chain_1)[2]
  
  ## params is the number of parameters (dimension of sampling space)
  ## that we have
  
  params <- dim(x$chain_1)[1]
  rhats <- numeric()
  for (i in 1:params) {
    
    ## I'll document more later
    ## this is just some variance / mean calculation 
    ## across and within the different chains
    
    
    mu_each_chain <- as.numeric(lapply (x, function(y) { mean(y[i,])})   )
    
    theta_2bar <- (1/m) * sum(mu_each_chain)
    
    B <- (n / (m-1)) * sum (  (mu_each_chain - theta_2bar)^2 )
    
    W <- sum(as.numeric(lapply(x, function(chain_matrix) {var(chain_matrix[i, ])}))) / m
    
    R2 <- (W * (1 - 1/n) + B / n) / W
    rhats <- c(rhats, sqrt(R2))
    
  }
  
  
  return(rhats)
  
  
  
  
}