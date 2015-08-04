context("Testing Hessian of Log Barrier")

test_that("Testing Hessian of Log Barrier", {
  
  ## Initialize toy A, b, x 
  
  A <- matrix(c(1,2,3,4), ncol = 2, nrow = 2, byrow = T)
  b <- c(10, 20)
  x <- c(1, 1)
  
  ## just to make sure x indeed satisfy Ax <= b
  
  expect_true(all( A %*% x <= b))
  
  ## compute the Hessian of the log barrier
  
  hess <- hessian_logbarrier(A = A, b = b, x = x)
  by.hand.hess <- t(A) %*% matrix(c(1/49, 0, 0, 1/169), ncol = 2, nrow = 2, byrow = T) %*% A
  
  ## they should be equal
  
  expect_true( all(hess - by.hand.hess <= 1e-10))
  
  ## hessian should be square, and symmetric, and have non-zero det
  
  expect_equal(ncol(hess), nrow(hess))
  expect_true(isSymmetric(hess))
  expect_true(det(hess) != 0)
  
  
})