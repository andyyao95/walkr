context("Testing Walkr Output")

## Tests that Walkr returns the desired 
## output

test_that("Testing Walkr Output", {
  
  ## Simple 3D simplex
  
  A1 <- matrix(c(1,0,1), ncol = 3)
  b1 <- 0.5
  
  ## Random 50D A 
  
  A2 <- matrix(sample(c(1,0), 50, replace = T), ncol = 50)
  b2 <- 0.6

  z1_har <- walkr(A = A1, b = b1, points = 50, method = "hit-and-run", chains = 5)
  z1_dikin <- walkr(A = A1, b = b1, points = 50, method = "dikin", chains = 5)
  
  z2_har <- walkr(A = A2, b = b2, points = 50, method = "hit-and-run", chains = 5)
  z2_dikin <- walkr(A = A2, b = b2, points = 50, method = "dikin", chains = 5)
  
  ## check that we're returning a list
  
  expect_equal(class(z1_har), "list")
  expect_equal(class(z1_dikin), "list")
  expect_equal(class(z2_har), "list")
  expect_equal(class(z2_dikin), "list")
  
  ## check that we have 5 chains in each
  
  expect_equal(length(z1_har), 5)
  expect_equal(length(z1_dikin), 5)
  expect_equal(length(z2_har), 5)
  expect_equal(length(z2_dikin), 5)
  
  ## check that each chain has a matrix
  
  for (i in 1:5) {
    
    expect_true(is.matrix(z1_har[[i]]))
    expect_true(is.matrix(z1_dikin[[i]]))
    expect_true(is.matrix(z2_har[[i]]))
    expect_true(is.matrix(z2_dikin[[i]]))
    
    
  }
  
  ## check that each chain has the right dimension
  
  for (i in 1:5) {
    
    expect_equal(dim(z1_har[[i]]), c(3, 50))
    expect_equal(dim(z1_dikin[[i]]), c(3, 50))
    expect_equal(dim(z2_har[[i]]), c(50, 50))
    expect_equal(dim(z2_dikin[[i]]), c(50, 50))
    
    
  }
  
  ## check that the names of each element of the list is correct
  
  expect_equal(names(z1_har), c("chain_1", "chain_2", "chain_3", "chain_4", "chain_5"))
  expect_equal(names(z1_dikin), c("chain_1", "chain_2", "chain_3", "chain_4", "chain_5"))
  expect_equal(names(z2_har), c("chain_1", "chain_2", "chain_3", "chain_4", "chain_5"))
  expect_equal(names(z2_dikin), c("chain_1", "chain_2", "chain_3", "chain_4", "chain_5"))
  
  ## check that we have no NA's in the output
  
  expect_true(!any(is.na(z1_har)))
  expect_true(!any(is.na(z1_dikin)))
  expect_true(!any(is.na(z2_har)))
  expect_true(!any(is.na(z2_dikin)))
  
})
