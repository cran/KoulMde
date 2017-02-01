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
List GLI(NumericVector E, int n, int SN, int EN) {
  List OutList(3);
  
  NumericVector NewE(n);
  
  for(int i=0; i<n; i++){
    NewE[i]=E[i];
  }
  
  int StartNum = SN;
  int EndNum = SN+(n-1);
  
  OutList[0] = NewE;
  OutList[1] = StartNum;
  OutList[2] = EndNum;
    
    
    
  return OutList  ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
