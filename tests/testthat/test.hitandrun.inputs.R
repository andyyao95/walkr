# ## In this file, we test different inputs to the hitandrun function
# 
# context("Testing hitandrun inputs")
# 
# test_that("sample points not a multiple of number of chains", {
#   A <- matrix(1, ncol = 3)
#   b <- 1
#   expect_error(hitandrun(A, b, n = 10, chains = 3))
# })
# 
# test_that("n must be positive integer", {
#   
#   ## n is not positive
#   
#   A <- matrix(1, ncol = 3)
#   b <- 1
#   expect_error(hitandrun(A, b, n = -1), "n must be a positive integer")
#   expect_error(hitandrun(A, b, n = 0),  "n must be a positive integer")
#   
#   ## n is not an integer
#   
#   expect_error(hitandrun(A, b, n = 1.5))
# })
# 
# test_that("if there are NAs in inputs of hitandrun",{
#   
#   ## NA's in b
#   
#   A = matrix(c(1, 1, 1), ncol = 3)
#   b = NA
#   expect_error(hitandrun(A, b, n = 1))
#   
#   ## NA's in A
#   
#   A = matrix(c(NA, NA, 1), ncol = 3)
#   b = 1
#   expect_error(hitandrun(A, b, n = 1))})
# 
# 
# test_that("Problems aren't overdetermined", {
#   A = matrix(rnorm(6), ncol = 2, nrow = 3)
#   b = c(1,2,3)
#   expect_error(hitandrun(A, b, n = 1))})  