#' Visualise Sampling
#' 
#' This function takes in a set of sample points and diagnoses them. The final
#' results are displayed in a Shiny app. The app contains the confidence
#' interval of each dimension's coordiates across the whole sets of points,
#' and draws a graph about the distribution of all points.
#' @param x is the set of points sampled 
#' 
#' @return a shiny interface that display the set of sampled points
#' @export
#' 

vis_sampling <- function(x, chains) {
      
      ## checking
      stopifnot(is.matrix(x))
      
      ## create local copy since we are performing side-effects
      
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