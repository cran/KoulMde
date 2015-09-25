#' Performs minimum distance estimation in linear regression model: Y=Xb + error
#'@param Y - Response variable in linear regression model
#'@param X - Explanatory variable in linear regression model
#'@return Returns betahat - Minimum distance estimator of b 
#'@examples
#'####################
#'n <- 10
#'p <- 3
#'X <- matrix(rnorm(n*p, 5,3), nrow=n, ncol=p)  #### Generate n-by-p design matrix X 
#'beta <- c(-2, 0.3, 1.5)                       #### Generate true beta = (-2, 0.3, 1.5)' 
#'eps <- rnorm(n, 0,1)                          #### Generate errors from N(0,1)
#'Y <- X%*%beta + eps
#'#####################
#'
#'betahat <- KoulLrMde(Y,X)                     ##### Obtain minimum distance estimator betahat 



#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@seealso KoulArMde()
#'@export
#'@importFrom "stats" "optim"




KoulLrMde <- function(Y, X){
  
  DimMat <- dim(X)
  LengY <- length(Y)  
  if (is.null(DimMat) == TRUE ){
    message("X should be a matrix")
    stop		
  }else{
    nXRow <- DimMat[1]
    nXCol <- DimMat[2]
    
    if (nXRow != LengY){
      message("Dimension of X does not match dimension of Y")
      stop
    }
  }
  
  BetaOLS <- solve(t(X)%*%X)%*% (t(X)%*%Y)
  Tmin <- optim(BetaOLS, TLRLoss(Y, X))
  return(Tmin$par)
}


TLRLoss <- function(Y, X){
  
  DimMat <- dim(X)
  nXRow <- DimMat[1]
  nXCol <- DimMat[2]
  
  A <- (t(X)%*%X)^(-1/2);
  
  Dual <- function(t){
    fval <- 0
    for (k in 1:nXCol){
      ak <- A[,k]
      for (i in 1:nXRow){
        xi <- t(X[i,])
        dik <- xi %*% ak
        ei <- Y[i] - xi %*% t
        for(j in i:nXRow){
          
          xj <- t(X[j,])
          djk <- xj %*% ak
          ej <- Y[j] - xj %*% t
          if (j==i){
            fval <- fval + djk^2* 2*abs(ej)
          }else{
            fval <- fval + dik*djk*Status(ei, ej)
          }
        }
        
      }
      
    }
    return(fval)
    
  }
  
  return(Dual)
}



Status <- function(a,b){
  if (a>=0 && a+b>=0){
    return( 2*(a+b))
  }else {
    return( -2*(a+b))
    
  }
  
}


#' Performs minimum distance estimation in autoregression model
#'@param X : vector of n observed value
#'@param AR_Order : oder of the autoregression model
#'@return returns minimum distance estimators of the parameter in the autoregression model
#'@examples
#'##### Generate stationary AR(2) process with 10 observations 
#'n <- 10
#'p <- 2
#'rho <- c(-0.2, 0.8)    ##### Generate true parameters rho = (-0.2, 0.8)'
#'eps <- rnorm(n, 0,1)   ##### Generate innovations from N(0,1)
#'X <- rep(0, times=n)
#'for (i in 1:n){
#'  tempCol <- rep(0, times=p)
#'  for (j in 1:p){
#'    if(i-j<=0){
#'      tempCol[j] <- 0
#'    }else{
#'      tempCol[j] <- X[i-j]
#'    }
#'  } 
#'X[i] <- t(tempCol)%*% rho + eps[i]
#'}
#'rhohat <- KoulArMde(X, p)

#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@importFrom "stats" "nlminb"
#'@export
#'@seealso KoulLrMde()














KoulArMde <- function(X, AR_Order){
  
  nLength <- length(X)
  
  if(nLength<=AR_Order){
    message("Length of vector X should be greater than AR_Order.")
    stop
  }
  
  Xres <- rep(0, times=(nLength-AR_Order))
  tempvec <- rep(0, times= AR_Order*(nLength-AR_Order) ) 
  Xexp <- matrix( tempvec, nrow = (nLength-AR_Order), ncol = AR_Order)
  
  for (i in 1:(nLength-AR_Order) ) {
    Xres[i] <- X[nLength - (i-1)]
    for (j in 1:AR_Order){
      Xexp[i,j] <- X[nLength-(i+j-1) ]
      
    }
  }               
  
  tempdet <- det(  t(Xexp) %*% Xexp )
  if (  tempdet < 0.01 ){
    rho0 <- 0.5*rep(1, times = AR_Order)  
  }else{
    rho0 <- solve(t(Xexp)%*%Xexp)%*% (t(Xexp)%*%Xres)
  }
  
  lbVec <- rep(-1, times=AR_Order)
  ubVec <- rep(1, times=AR_Order)
  
  Tmin <- nlminb(rho0, TARLoss(X, AR_Order), lower=lbVec, upper=ubVec)
  
  return(Tmin$par)
}


TARLoss <- function(X, AR_Order){
  
  nLength <- length(X)
  
  Dual <- function(r){
    fval <- 0
    for (k in 1:AR_Order){
      
      for (i in 1:nLength){
        if(i<=k){
          dik <- 0
        }else{
          dik <- X[i-k]/sqrt(nLength)
        }
        
        tempColi <- rep(0, times=AR_Order)
        for(m in 1:AR_Order){
          if(i-m<=0){
            tempColi[m] <- 0
          }else{
            tempColi[m] <- X[i-m]
          }	
          
        }
        ei <- X[i] - t(r)%*%tempColi 
        
        for(j in i:nLength){
          if(j<=k){
            djk <- 0
          }else{
            djk <- X[j-k]/sqrt(nLength)
          }
          tempColj <- rep(0, times=AR_Order)
          
          for(m in 1:AR_Order){
            if(j-m<=0){
              tempColj[m] <- 0
            }else{
              tempColj[m] <- X[j-m]
            }	
            
          }
          ej <- X[j] - t(r)%*%tempColj
          
          if (j==i){
            fval <- fval + djk^2* 2*abs(ej)
          }else{
            fval <- fval + dik*djk*Status(ei, ej)
          }
          
        }
        
      }
      
    }
    return(fval)
    
  }
  
  return(Dual)
}


