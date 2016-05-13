#' Performs two-stage minimum distance estimation in linear regression model with autoregressive error.
#'@param Y - Vecotor of response variable in linear regression model.
#'@param X - Design matrix of explanatory variable in linear regression model.
#'@param D - Weight Matrix. Dimesion of D should match that of X. Default value is X*A where A=(X'*X)^(-1/2).
#'@param RegIntMeasure - Measure used in integration for the estimation of regression coefficients. It should be either Lebesgue or degenerate. 
#'@param AR_Order - Oder of the autoregressive error.
#'@param ArIntMeasure - Measure used in the first stage for estimation of autoregressive coefficient of error.
#'@return MDE1stage - The list of the first stage minimum distance estimation result. It contains betahat1stage, residual1stage, and rho1stage.    
#'@return betahat1stage - The first stage minimum distance estimator of regression coefficient.
#'@return residual1stage - Residuals after the first stage minimum distance estimation. 
#'@return rho1stage - The first stage minimum distance estimator of autoregressive coefficient.
#'@return MDE2stage - The list of the second stage minimum distance estimation result. It contains betahat2stage, residual2stage, and rho2stage.
#'@return betahat2stage - The second stage minimum distance estimator of regression coefficient.
#'@return residual2stage - Residuals after the second stage minimum distance estimation. 
#'@return rho2stage - The second stage minimum distance estimator of autoregressive coefficient.
#'@examples
#'####################
#'n <- 10
#'p <- 3
#'X <- matrix(rnorm(n*p, 5,3), nrow=n, ncol=p)   #### Generate n-by-p design matrix X 
#'beta <- c(-2, 0.3, 1.5)                        #### Generate true beta = (-2, 0.3, 1.5)' 
#'rho  <- 0.4                                    #### True rho = 0.4  
#'eps <- vector(length=n)                        
#'xi <- rnorm(n, 0,1)                            #### Generate innovation from N(0,1)
#'                                               #### Generate autoregressive process or order 1
#'for(i in 1:n){       
#'  if(i==1){eps[i] <- xi[i]}
#'  else{eps[i] <- rho*eps[i-1] + xi[i]}
#'} 
#'Y <- X%*%beta + eps
#'#####################
#'D <- "default"                                  #### Use the default weight matrix 
#'MDEResult <- Koul2StageMDE(Y,X, "default", "Lebesgue", 1, "Lebesgue")
#'MDE1stageResult <- MDEResult[[1]]
#'MDE2stageResult <- MDEResult[[1]]
#'
#'beta1 <- MDE1stageResult$betahat1stage
#'residual1 <- MDE1stageResult$residual1stage
#'rho1 <- MDE1stageResult$rhohat1stage
#'
#'beta2 <- MDE2stageResult$betahat1stage
#'residual2 <- MDE1stageResult$residual2stage
#'rho2 <- MDE2stageResult$rhohat1stage



#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@seealso KoulArMde() and KoulLrMde() 
#'@export
#'@importFrom "stats" "optim"
#'@importFrom "stats" "nlm"

Koul2StageMDE <- function(Y,X,D,RegIntMeasure, AR_Order, ArIntMeasure){
  
  
  DimMat <- dim(X)
  n <- DimMat[1]
  p <- DimMat[2]
  
  MDE1Result <- KoulLrMde(Y,X, D, RegIntMeasure)
  betahat1stage <- MDE1Result$betahat
  residual1stage <- MDE1Result$residual
  rho1stage <- KoulArMde(residual1stage, AR_Order, ArIntMeasure)
  
  MDE1stage <- list(betahat1stage, residual1stage, rho1stage)
  
  ###########################   2 stage MDE
  Ytilde <- vector(length=(n-1))
  Xtilde <- matrix(rep(0,times=(n-1)*p), nrow=(n-1), ncol=p )
  
  for(j in 1:(n-1)){
    Xtilde[j, ] <- X[(j+1), ] - rho1stage*X[j, ]
    Ytilde[j] <- Y[j+1] - rho1stage*Y[j]
  }
  MDE2Result <- KoulLrMde(Ytilde, Xtilde, D, RegIntMeasure)
  betahat2stage <- MDE2Result$betahat
  residual2stage <- Y-X%*%betahat2stage
  rho2stage <- KoulArMde(residual2stage, AR_Order, ArIntMeasure)
  
  MDE2stage <- list(betahat2stage, residual2stage, rho2stage)
  
  ResultVal <- list(MDE1stage, MDE2stage)
  return(ResultVal)
  
}




#' Performs minimum distance estimation in linear regression model: Y=Xb + error.
#'@param Y - Vecotor of response variable in linear regression model.
#'@param X - Design matrix of explanatory variable in linear regression model.
#'@param D - Weight Matrix. Dimesion of D should match that of X. Default value is X*A where A=(X'*X)^(-1/2).
#'@param IntMeasure - Measure used in integration. It should be either Lebesgue or degenerate. 
#'@return betahat - Minimum distance estimator of b. 
#'@return residual - Residuals after minimum distance estimation. 
#'@examples
#'####################
#'n <- 10
#'p <- 3
#'X <- matrix(rnorm(n*p, 5,3), nrow=n, ncol=p)   #### Generate n-by-p design matrix X 
#'beta <- c(-2, 0.3, 1.5)                        #### Generate true beta = (-2, 0.3, 1.5)' 
#'eps <- rnorm(n, 0,1)                           #### Generate errors from N(0,1)
#'Y <- X%*%beta + eps
#'#####################
#'D = "default"                                  #### Use the default weight matrix 
#'KM <- KoulLrMde(Y,X,D,"Lebesgue")              ##### Use Lebesgue measure for integration
#'betahat <- KoulLrMde(Y,X,D,"Lebesgue")$betahat ##### Obtain minimum distance estimator 
#'resid <- KoulLrMde(Y,X,D,"Lebesgue")$residual  ##### Obtain residual 


#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@seealso KoulArMde() and Koul2StageMDE()
#'@export
#'@importFrom "stats" "optim"
#'@importFrom "stats" "nlm"



KoulLrMde <- function(Y, X, D, IntMeasure){
  
  if(nargs() != 4){
    message("Number of arguments should be four.")
    stop()
  }
  Hx = IntMeasure
  if( (Hx != "Lebesgue") && (Hx != "Degenerate")){
    message("Intergrating measure Hx should be either Lebesgue or Degenerate.")
    stop()
  }

  if (is.vector(X) == TRUE ){
    
    nXRow <- length(X)
    nXCol <- 1
    LengY <- length(Y)

    if (nXRow != LengY){
      message("Dimension of X does not match dimension of Y.")
      stop()
    }
    
    if(is.vector(D) == TRUE){
      nDRow <-  length(D)
      nDCol <- 1
      
    }else{
      message("When X is a vector, D should be a vector too.")
      stop()
    }
    
    if(nDRow != nXRow){
      str= paste("D should be ", nXRow, "-by-1 vector.")
      message(str)
      stop()
    }
  
  }else if(is.matrix(X) == TRUE){
    DimMat <- dim(X)
    
    LengY <- length(Y)  
    
    nXRow <- DimMat[1]
    nXCol <- DimMat[2]
    
    if(is.matrix(D) == TRUE){
      DDimMat <- dim(D)
    }else if(D == "default"){
      tempA <- (t(X)%*%X)
      A <- sqrtmat(tempA, -0.5)
      D <- X%*%A
      
    }else{
      message("D should be a matrix.")
      stop()
    } 
    
    DDimMat <- dim(D)
    nDRow <- DDimMat[1]
    nDCol <- DDimMat[2]
    
    if (nXRow != LengY){
      message("Dimension of X does not match dimension of Y.")
      stop()
    }
    
    if( (nXRow != nDRow) || ((nXCol != nDCol)) ) {
      message("Dimesion of D should match dimension of X.")
      stop()
    }
    
  }else{
    message("X is not a valid design matrix.")
    stop()
    
  }

  
  if (nXCol == 1){
    BetaOLS = solve(t(X)%*%X)%*% (t(X)%*%Y)
    Tmin <- nlm(TLRLoss(Y, X, D, Hx), BetaOLS)
    bhat <- Tmin$estimate
  }else{
    BetaOLS <- solve(t(X)%*%X)%*% (t(X)%*%Y)
    Tmin <- optim(BetaOLS, TLRLoss(Y, X, D, Hx))
    bhat <- Tmin$par
  }
  
  if (is.vector(X) == TRUE ){
    res <- Y - bhat*X
  }else{
    res <- Y - X %*% bhat  
  }
  
  lst = list(betahat=bhat, residual = res)
  return(lst)
  
}


sqrtmat<-function(MAT, EXP, tol=NULL){
  MAT <- as.matrix(MAT)
  matdim <- dim(MAT)
  if(is.null(tol)){
    tol=min(1e-7, .Machine$double.eps*max(matdim)*max(MAT))
  }
  if(matdim[1]>=matdim[2]){ 
    svd1 <- svd(MAT)
    keep <- which(svd1$d > tol)
    res <- t(svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep]))
  }
  if(matdim[1]<matdim[2]){ 
    svd1 <- svd(t(MAT))
    keep <- which(svd1$d > tol)
    res <- svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep])
  }
  return(res)
}




TLRLoss <- function(Y, X, D, Hx){
  
  if (is.vector(X) == TRUE ){
    nXRow <- length(X)
    nXCol <- 1

  }else {
    DimMat <- dim(X)
    LengY <- length(Y)  
    
    nXRow <- DimMat[1]
    nXCol <- DimMat[2]

  }

  
  Dual <- function(t){
    fval <- 0
    for (k in 1:nXCol){
      if(Hx == "Lebesgue"){
        
        for (i in 1:nXRow){
          if(nXCol == 1){xi <- X[i]} else {xi <- t(X[i,])}  
          
          dik <- D[i,k]
          ei <- Y[i] - xi %*% t
          for(j in i:nXRow){
            if(nXCol == 1){xj <- X[j]} else {xj <- t(X[j,])}
            
            djk <- D[j,k]
            ej <- Y[j] - xj %*% t
            
            fval <- fval + 2*dik*djk* ( max(ei, -ej)+max(-ei, ej) - max(ei, ej)-max(-ei, -ej))
            
          }
        }
      }else{
        tempVal <- 0
        
        for(i in 1:nXRow){
          dik <- D[i,k]
          if(nXCol == 1){xi <- X[i]} else {xi <- t(X[i,])}
          ei <- Y[i] - xi %*% t
          if(ei<0){
            sgn <- -1
          }else if(ei == 0){
            sgn <- 0
          }else{
            sgn <- 1
          }
          tempVal <- tempVal + dik*sgn
        }
        fval <- fval + tempVal^2
      }
    }
    return(fval)
    
  }
  
  return(Dual)
}



#' Performs minimum distance estimation in autoregression model.
#'@param X : vector of n observed value.
#'@param AR_Order : oder of the autoregression model.
#'@param IntMeasure - Measure used in integration. It should be either Lebesgue or degenerate.
#'@return returns minimum distance estimators of the parameter in the autoregression model.
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
#'rhohat <- KoulArMde(X, p, "Lebesgue")

#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@importFrom "stats" "nlminb"
#'@export
#'@seealso KoulLrMde() and Koul2StageMDE()





KoulArMde <- function(X, AR_Order, IntMeasure){
  
  Hx = IntMeasure
  if( (Hx != "Lebesgue") && (Hx != "Degenerate")){
    message("Intergrating measure Hx should be either Lebesgue or Degenerate.")
    stop()
  }
  
  nLength <- length(X)
  
  if(nLength<=AR_Order){
    message("Length of vector X should be greater than AR_Order.")
    stop()
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
  
  Tmin <- nlminb(rho0, TARLoss(X, AR_Order, Hx), lower=lbVec, upper=ubVec)
  
  return(Tmin$par)
}


TARLoss <- function(X, AR_Order, Hx){
  
  nLength <- length(X)
  
  Dual <- function(r){
    fval <- 0
    for (k in 1:AR_Order){
      
      if(Hx == "Lebesgue"){
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
            
            fval <- fval + 2*dik*djk* ( max(ei, -ej)+max(-ei, ej) - max(ei, ej)-max(-ei, -ej))
          }
          
        }
        
        
      }else{
        tempVal<-0
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
          
          if(ei<0){
            sgn <- -1
          }else if(ei == 0){
            sgn <- 0
          }else{
            sgn <- 1
          }
          tempVal <- tempVal+dik*sgn
        }
        fval <- fval+tempVal^2
      }
      

    }
    return(fval)
    
  }
  
  return(Dual)
}


