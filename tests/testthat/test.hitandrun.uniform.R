context("Testing hit-and-run uniformity")

test_that("Testing hit-and-run uniformity", {
  
  ## Simplest possible case 
  
  A <- matrix(1, ncol = 2)
  b <- 1
  
  ## In this case, we are sampling from the 2D simplex
  
  ## Draw a sample of 1,000. If the sample is truly random,
  ## then each variable should have a 50/50 chance of being above 0.5. It 
  ## follows that the number of times, out of a thousand, that it is above 0.5 
  ## is provided by a binomial distrbution. So, we check that number of times 
  ## above 0.5, for both variables, is less than the 99th percentile of the
  ## binomial distribution.
  
  z <- walkr(A = A, b = b, n = 1000, method = "hit-and-run")
  expect_true(all(apply(z, 1, function(x) length(which(x> .5)) < 
                          qbinom(.99, 1000, .5))))
  
  ## also should be greater than the 1st percentile of the binomial
  
  expect_true(all(apply(z, 1, function(x) length(which(x> .5)) > 
                          qbinom(.01, 1000, .5))))
  
 
  ### What about a 5D Simplex
  
  A <- matrix(1, ncol = 5)
  b <- 1 
  
  z <- walkr(A = A, b = b, n = 5000, method = "hit-and-run")
  
  ## should expect that the sum of x_1, x_2 be roughly the same as x_4, x_5
  
  sum1 <- sum(z[1,]) + sum(z[2,])
  sum2 <- sum(z[4,]) + sum(z[5,])
  
  ## expected value of their difference is zero
  ## normal approximation applicable here
  ## 99% confidence interval
  
  conf_interval.99 <- qnorm(p = c(0.01, 0.99), mean = 0, sd = sqrt(5000*0.5*0.5))
  
  expect_true(sum1-sum2 >= conf_interval.95[1])
  expect_true(sum1-sum2 <= conf_interval.95[2])
   
})