// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#ifdef __GLIBC__
# define _POSIX_C_SOURCE 200809L
#endif
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;
using Eigen::MatrixXd;

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

// ONLY WRITE C CODE / RCPP FUNCTIONS IN THIS FILE BELOW
// USE THIS RCPP MAGIC TO DEAL WITH EVERYTHING ELSE =====> Rcpp::compileAttributes()

// 1. The transpose of a matrix
//
//' Fast Matrix Transpose
//' 
//' Computes using RcppEigen transposed of A
//' 
//' @param A is the matrix being transposed
//' 
//' @return transpose of A
//' @examples
//' \dontrun{
//' rcppeigen_ftrans(A)
//' }
//' 
// [[Rcpp::export]]

Eigen::MatrixXd rcppeigen_ftrans(const Eigen::MatrixXd & A) {
  Eigen::MatrixXd m = A.transpose();
  return m;
}


// 2. solve(A)

//' Fast Matrix Inverse
//' 
//' Computes using RcppEigen the inverse of A
//' 
//' @param A is the matrix being inverted
//' 
//' @return inverse of A
//' 
//' @examples
//' \dontrun{
//' rcppeigen_fsolve(A)
//' }
//' 
//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_fsolve(const Eigen::Map<Eigen::MatrixXd> & A){
  Eigen::MatrixXd Ainv = A.inverse();
  return Ainv;
}

// 3. det(A)

//' Fast Matrix Determinant
//' 
//' Computes using RcppEigen the determinant of A
//' 
//' @param A is the matrix whose determinant calculated 
//' 
//' @return determinant of A
//' 
//' @examples
//' \dontrun{
//' rcppeigen_fdet(A)
//' }
//' 
//[[Rcpp::export]]
double rcppeigen_fdet(const Eigen::Map<Eigen::MatrixXd> & A){
 return A.determinant();
}

// 4. A %*% B

//' Fast Matrix Product
//' 
//' Computes using RcppEigen the product of A and B
//' 
//' @param A is the first parameter in A times B
//' @param B is the second parameter in A times B
//' 
//' @return matrix product A times B
//' 
//' @examples
//' \dontrun{
//' rcppeigen_fprod(A, B)
//' }
//' 
//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_fprod(const Eigen::Map<Eigen::MatrixXd> & A, 
                       const Eigen::Map<Eigen::MatrixXd> & B){
 return A * B;
}

// 5. t(A) %*% B

//' Fast Matrix Cross-Product
//' 
//' Computes using RcppEigen the product of t(A) and B
//' 
//' @param A is the first parameter in t(A) times B
//' @param B is the second parameter in t(A) times B
//' 
//' @return matrix cross-product t(A) times B
//' 
//' @examples
//' \dontrun{
//' rcppeigen_fcrossprod(A, B)
//' }
//' 
//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_fcrossprod(const Eigen::Map<Eigen::MatrixXd> & A, 
                       const Eigen::Map<Eigen::MatrixXd> & B){
 return A.transpose() * B;
}

// 6. A %*% t(B)

//' Fast Matrix T-Cross-Product
//' 
//' Computes using RcppEigen the product of A and t(B)
//' 
//' @param A is the first parameter in A times t(B)
//' @param B is the second parameter in A times t(B)
//' 
//' @return matrix tcross-product A times t(B)
//' 
//' @examples
//' \dontrun{
//' rcppeigen_ftcrossprod(A, B)
//' }
//' 
//[[Rcpp::export]]
Eigen::MatrixXd rcppeigen_ftcrossprod(const Eigen::Map<Eigen::MatrixXd> & A, 
                       const Eigen::Map<Eigen::MatrixXd> & B){
 return A * B.transpose();
}