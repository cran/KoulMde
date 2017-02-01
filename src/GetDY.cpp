#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix DY(NumericVector Y, NumericMatrix D){
  int nrow = D.nrow();
  int ncol = D.ncol();
  double SumD = 0;
  
  NumericMatrix OutMatrix(3*nrow, nrow);
  
  for(int i=0; i<nrow; i++){
    for(int j=0; j<nrow; j++){
      SumD=0;
      for(int k=0; k<ncol; k++){
        SumD += D(i,k)*D(j,k);
      }
      OutMatrix(i,j) = SumD;
      OutMatrix((i+nrow),j)=Y[i]+Y[j];
      OutMatrix((i+2*nrow),j)=Y[i]-Y[j];
    }
  }
  return OutMatrix ;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

//
//timesTwo(42)
//
