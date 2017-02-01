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
NumericMatrix Xpm(NumericMatrix X){
  //NumericVector Xpm(NumericMatrix X){  
  
  int nrow = X.nrow();
  int ncol = X.ncol();
  
  NumericMatrix outMatrix(2*nrow*nrow, ncol);

  NumericVector Xi(ncol);
  NumericVector Xj(ncol);

  
  for(int i=0; i<nrow; i++){
    

    for(int j=0; j<nrow; j++){
      
      for(int k=0; k<ncol; k++){
        outMatrix( (i*nrow)+j, k) = X(i,k) + X(j,k);
        outMatrix( (i*nrow)+j+(nrow*nrow), k) = X(i,k) - X(j,k);
      }
      
    }
    
    
  }
  
  
  return outMatrix;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
