context("Testing dikin uniformity")

test_that("Testing dikin uniformity", {
  
  ## The Simple 3D simplex
  
  A <- matrix(1, ncol = 3)
  b <- 1
  
  ## We should expect to see that all coordinates have
  ## the same value. Note that the expected value of each of the coordinates 
  ## is 1/3 
  
  ## construct confidence interval
  
  conf <- qnorm(p = c(0.01, 0.99), mean = 1/3, sd = sqrt(0.5*0.5/1000))
  
  z <- walkr(A = A, b = b, n = 1000, method = "dikin")
  
  ## all should fall within confidence interval
  
  expect_true(mean(z[1,]) <= conf[2])
  expect_true(mean(z[1,]) >= conf[1])
  expect_true(mean(z[2,]) <= conf[2])
  expect_true(mean(z[2,]) >= conf[1])
  expect_true(mean(z[3,]) <= conf[2])
  expect_true(mean(z[3,]) >= conf[1])
  
  
  
})