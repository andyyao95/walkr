#' Visualise Sampling
#' 
#' This function takes in a set of sample points and diagnoses them using
#' the shinyStan interface. The app contains the confidence
#' interval of each dimension's coordiates across the whole sets of points,
#' Gelman-Rubin statistics, trace-plots, and other diagnostic tools for
#' examining convergence.
#' 
#' @param x is the set of points sampled, with its columns as the sampled points.
#'        If multiple chains are present, then the columns should be ordered
#'        such that each chain follow each other. 
#' @param chains is the number of chains that x contains
#' 
#' @return a shiny interface that display the diagnostics of the MCMC random walk
#' @export
#' @importFrom shinyStan launch_shinystan as.shinystan
#' 

vis_sampling <- function(x, chains) {
      
      this.df <- x
      
      ## iterations is the number of points in hitandrun
      ## parameters is the dimension of the sample space
      
      iterations <- ncol(this.df) / chains
      parameters <- nrow(this.df)
      
      ## we have to transpose here in order to make the dimensions fit 
      
      this.df <- t(this.df)
      
      ## remake the data set into 3 dimensional array
      
      dim(this.df) <- c(iterations, chains, parameters)
      
      ## using the stan package
      
      shinyStan::launch_shinystan(shinyStan::as.shinystan(this.df))
      
}