context("Testing Dikin Ellipsoid Function")

test_that("Testing Dikin Ellipsoid Function", {
  
  ## Initialize toy A, b, x 
  
  A <- matrix(c(1,2,3,4), ncol = 2, nrow = 2, byrow = T)
  b <- c(10, 20)
  x <- c(1, 1)
  
  ## should be TRUE by hand
  expect_true(dikin_ellipsoid(A = A, b = b, x0 = x, r = 1, z = c(1.5,3)))
  
  ## should be FALSE by hand
  expect_true(!dikin_ellipsoid(A = A, b = b, x0 = x, r = 1, z = c(2,3)))

  
  
})