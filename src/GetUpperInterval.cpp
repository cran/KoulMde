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
List GUI(NumericVector E, int n, int SN, int EN) {
  
  int len = E.length();
  List OutList(3);
  
  int NewElen = (len-n);
  
  NumericVector NewE(NewElen);
  
  for(int i=0; i<NewElen; i++){
    NewE[i]=E[n+i];
  }
  
  int StartNum = EN-NewElen+1;
  int EndNum = EN;
  
  OutList[0] = NewE;
  OutList[1] = StartNum;
  OutList[2] = EndNum;
  
  
  
  return OutList  ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
