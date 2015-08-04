context("Testing Dikin Walk")

test_that("Testing Dikin Walk", {
  
  ## 3D simplex, initialize Aw = b  
  n <- 3
  A <- matrix(rep(1,n), ncol = n, nrow = 1)
  b <- 1 
  
  
  ## Find the basis representation
  ## performing affine transformation into 
  ## new A and new b (Ax <= b)
  
  z <- complete_solution(A,b, randomize = T)
  particular  <- z$particular
  homogeneous <- z$homogeneous
  new_A <- -homogeneous
  new_b <- particular
  
  
  ## find a starting point in the polytope
  
  x0 <- start_point(A = new_A, b = new_b, n = 1, average = 30)
  my.center <- as.vector(x0)
  
  ## run the algorithm
  
  z <- dikin_walk(A = new_A, b = new_b, n = 50, r = 1, x0 = my.center)
  answer <- apply( z, 2, function(x) { homogeneous %*% x + particular  })
  
  ## check results, indeed on the simplex
  expect_true(all(apply(answer, 2, sum) - 1 < 1e-10))
  
  ## non-negativity
  expect_true(length(which(answer < 0)) == 0)
  
  ######## HIGHER DIM WITH CONSTRAINTS ########
  
  ## 20D simplex, intersecting 3 hyperplanes, initialize Aw = b  
  n <- 20
  A <- matrix(rep(1,n), ncol = n, nrow = 1)
  A <- rbind(A, sample(c(1,0), n, replace = T))
  A <- rbind(A, sample(c(1,0), n, replace = T))
  A <- rbind(A, sample(c(1,0), n, replace = T))
  b <- c(1, 0.7, 0.2, 0.05)
  
  ## Find the basis representation
  ## performing affine transformation into 
  ## new A and new b (Ax <= b)
  
  z <- complete_solution(A,b, randomize = T)
  particular  <- z$particular
  homogeneous <- z$homogeneous
  new_A <- -homogeneous
  new_b <- particular
  
  
  ## find a starting point in the polytope
  
  x0 <- start_point(A = new_A, b = new_b, n = 1, average = 30)
  my.center <- as.vector(x0)
  
  ## run the algorithm
  
  z <- dikin_walk(A = new_A, b = new_b, n = 50, r = 1, x0 = my.center)
  answer <- apply( z, 2, function(x) { homogeneous %*% x + particular  })
  
  ## check results, indeed on the simplex
  expect_true(all(apply(answer, 2, sum) - 1 < 1e-10))
  
  ## non-negativity
  expect_true(length(which(answer < 0)) == 0)
 
  ## check constraints satisfied
  
  for (vect in 1:50) {
    
    expect_true(all((A %*% answer[,vect] - b) <= 1e-10 ))
    
  }
  
})
