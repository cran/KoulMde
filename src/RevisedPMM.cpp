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
NumericMatrix PMM(int n, NumericMatrix DYM, NumericMatrix XpmM, int l, int p, NumericVector bVec) {
  
  NumericMatrix outMat((2*n*n+1), 3);
  int PMLen = 0;
  int MMLen = 0;
  double tempSum = 0;
  
  double Xp;
  double Xm;
  

  for(int i=0; i<n; i++){

    for(int j=0; j<n; j++){
      
      Xp = XpmM((i*n)+j, l);
      
      if((Xp> 0) || (Xp<0)){
        
        if(Xp<0){
          outMat(PMLen,0) = 0-Xp;
        }else{
          outMat(PMLen,0) = Xp;
        }
        outMat(PMLen,1) = DYM(i,j);
        
        tempSum = 0;
        for(int k=0; k<p; k++){
          tempSum += bVec[k]* XpmM((i*n)+j, k);
        }
        outMat(PMLen,2) = (DYM(n+i,j) - tempSum + Xp*bVec[l])/Xp;  
        
        PMLen += 1;
        
      }
      
      Xm = XpmM( (i*n)+j+(n*n), l);
      if( (Xm> 0) || (Xm<0)  ){
        
        if(Xm<0){
          outMat((n*n)+MMLen,0) = 0-Xm;
        }else{
          outMat((n*n) + MMLen,0) = Xm;
        }
        outMat((n*n) + MMLen,1) = DYM(i,j);
        
        tempSum = 0;
        for(int k=0; k<p; k++){
          tempSum += bVec[k]*XpmM( (i*n)+j+(n*n), k);
        }
        outMat((n*n) + MMLen,2) = (DYM(2*n+i,j) - tempSum + Xm*bVec[l])/Xm;  
        MMLen += 1;
      }
      
      
      
      
    }
  }
  outMat((2*n*n), 0) = PMLen;
  outMat((2*n*n), 1) = MMLen;
  
  return outMat;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

