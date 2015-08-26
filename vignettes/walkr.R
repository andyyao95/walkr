## ----example1, eval = TRUE, cache = TRUE---------------------------------
A <- matrix(1, ncol = 3)
b <- 1
set.seed(314)
sampled_points <- walkr(A = A, b = b, points = 1000, 
                        method = "hit-and-run", chains = 5, ret.format = "matrix")

## ----example3, eval = TRUE, cache = TRUE---------------------------------
A <- matrix(sample(c(0,1,2), 40, replace = TRUE), ncol = 20)
b <- c(0.5, 0.3)
sampled_points <- walkr(A = A, b = b, points = 10000, chains = 5, 
                        method = "hit-and-run", ret.format = "list")


## ----example4, eval = TRUE, cache = TRUE---------------------------------
sampled_points <- walkr(A = A, b = b, points = 1000, chains = 5, thin = 500, 
                        method = "hit-and-run", ret.format = "list")         

## ----example5, eval = TRUE, cache = TRUE---------------------------------
sampled_points <- walkr(A = A, b = b, points = 1000, chains = 5, thin = 10,  
                        method = "dikin", ret.format = "list")

