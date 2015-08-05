// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// 1. The transpose of a matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd rcppeigen_ftrans(const Eigen::MatrixXd & A) {
  Eigen::MatrixXd m = A.transpose();
  return m;
}


// 2. solve(A)

//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_fsolve(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::MatrixXd Ainv = A.inverse();
  return Ainv;
}

// 3. det(A)

//[[Rcpp::export]]
double rcppeigen_fdet(const Eigen::Map<Eigen::MatrixXd> & A){
 return A.determinant();
}

// 4. A %*% B

//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_fprod(const Eigen::Map<Eigen::MatrixXd> & A, 
                       const Eigen::Map<Eigen::MatrixXd> & B){
 return A * B;
}

// 5. t(A) %*% B

//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_fcrossprod(const Eigen::Map<Eigen::MatrixXd> & A, 
                       const Eigen::Map<Eigen::MatrixXd> & B){
 return A.transpose() * B;
}

// 6. A %*% t(B)

//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_ftcrossprod(const Eigen::Map<Eigen::MatrixXd> & A, 
                       const Eigen::Map<Eigen::MatrixXd> & B){
 return A * B.transpose();
}
