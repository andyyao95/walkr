# walkr
**walkr** uses random walks to sample points from the intersection of 
the \eqn{N} simplex with $M$ hyperplanes. Mathematically, the sampling space is all vectors $x$ 
that satisfy $Ax=b$, $\sum x = 1$, and $x_i \geq 0$. The sampling algorithms implemented 
are hit-and-run and Dikin walk, both of which are MCMC (Monte-Carlo Markov Chain) random 
walks. **walkr** also provide tools to examine the convergence
properties of the random walks. 

# Getting Started

<!--   * Install from CRAN:

  `install.packages("walkr")`
-->
* Install from GitHub:  

  `devtools::install_github("andyyao95/walkr")`  

# Sampling Points  

  `library(walkr)`  
  `A <- matrix(1, ncol = 3)`    
  `b <- 1`    
  `sampled_points <- walkr(A = A, b = b, points = 1000)`    
  
# Visualizing the Sampled Points  

  `explore_walkr(sampled_points)`