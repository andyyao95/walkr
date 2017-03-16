// We are incredibly grateful to the Rcpp team for providing this code below 
// One can see the Rcpp github page for the thread discussing this issue with the CRAN Note

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_walkr(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
