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

double CLoss(NumericVector Y, NumericMatrix X, NumericMatrix Dstar, NumericVector beta){
  double total = 0;
  
  int nrow = X.nrow();
  int ncol = X.ncol();
  
  double tempxbi = 0;
  double tempxbj = 0;
  double absipj = 0;
  double absimj = 0;
  double ei = 0;
  double ej = 0;
  
  for(int i=0; i<nrow; i++){
    
    tempxbi=0;
    for(int k=0; k<ncol; k++){
      tempxbi += X(i,k)*beta[k];
    }
    
    ei = Y[i]- tempxbi;
      
    for(int j=i; j<nrow; j++){
      tempxbj=0;
      for(int h=0; h<ncol; h++){
        tempxbj += X(j,h)*beta[h];
      }
      
      ej = Y[j]- tempxbj;
      
      if((ei+ej)<0){
        absipj = 0-(ei+ej);
      }else{
        absipj = (ei+ej);
      }
      
      if((ei-ej)<0){
        absimj = 0-(ei-ej);
      }else{
        absimj = (ei-ej);
      }
      
      if( j== i){
        total += Dstar(i,j)* (absipj-absimj);
      }else{
        total += 2*Dstar(i,j)* (absipj-absimj);
      }
        
    }
    
  }
  return total;
  
}





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

