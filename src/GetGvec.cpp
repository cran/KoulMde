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
NumericVector GVec(NumericMatrix PMat) {
  
  int pmlen=0;
  pmlen = PMat.nrow();
  int pmlenp1 = pmlen+1;
  
  NumericVector OutVec(pmlenp1);
  
  double temp=0;
  
  for(int i=0; i<pmlenp1; i++){
    if(i==0){
      temp = 0;
      for(int j=0; j<pmlen; j++){
        temp += PMat(j, 0)*PMat(j, 1);
      }
      OutVec[i] = 0-( temp );
    }else{
      OutVec[i] = OutVec[i-1]+2*PMat((i-1), 0)*PMat((i-1), 1);  
    }
    
  }
  
  return OutVec ;
}




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
