#'Two-stage minimum distance estimation in linear regression model with autoregressive error.
#' 
#'Estimates both regression and autoregressive coefficients in the model \eqn{Y=X\beta + \epsilon} where \eqn{\epsilon} is autoregressive process of known order \code{q}
#'@param Y - Vector of response variables in linear regression model.
#'@param X - Design matrix of explanatory variables in linear regression model.
#'@param D - Weight matrix. Dimension of D should match that of X. "default" uses XA where A=(X'X)^(-1/2).
#'@param RegIntMeasure - Symmetric and \eqn{\sigma}-finite measure used for estimating \eqn{\beta}.  
#'@param AR_Order - Order of the autoregressive error.
#'@param ArIntMeasure - Symmetric and \eqn{\sigma}-finite measure used for estimating autoregressive coefficients of the error.
#'@return MDE1stage - The list of the first stage minimum distance estimation result. It contains betahat1stage, residual1stage, and rhohat1stage.    
#'\itemize{
#'  \item betahat1stage - The first stage minimum distance estimators of regression coefficients.
#'  \item residual1stage - Residuals after the first stage minimum distance estimation. 
#'  \item rhohat1stage - The first stage minimum distance estimators of autoregressive coefficients of the error.  
#'}
#'@return MDE2stage - The list of the second stage minimum distance estimation result. It contains betahat2stage, residual2stage, and rhohat2stage.
#'\itemize{
#'  \item betahat2stage - The second stage minimum distance estimators of regression coefficients.
#'  \item residual2stage - Residuals after the second stage minimum distance estimation. 
#'  \item rhohat2stage - The second stage minimum distance estimators of autoregressive coefficients of the error.
#'}
#'@examples
#'####################
#'n <- 10
#'p <- 2
#'X <- matrix(runif(n*p, 0,20), nrow=n, ncol=p)   #### Generate n-by-p design matrix X 
#'beta <- c(-2, 1.5)                        #### Generate true beta = (-2, 1.5)' 
#'
#'q <- 1
#'rho <- 0.8                            ##### Generate true parameters rho = 0.8
#'eps <- vector(length=n)                        
#'xi <- rnorm(n, 0,1)                            #### Generate innovation from N(0,1)
#'                                               #### Generate autoregressive process of order q=1
#'for (i in 1:n){
#'  tempCol <- rep(0, times=q)
#'  for (j in 1:q){
#'    if(i-j<=0){
#'      tempCol[j] <- 0
#'    }else{
#'      tempCol[j] <- eps[i-j]
#'    }
#'  } 
#'  eps[i] <- t(tempCol)%*% rho + xi[i]
#'}
#'
#'Y <- X%*%beta + eps
#'#####################
#'D <- "default"                                  #### Use the default weight matrix 
#'
#'Lx <- function(x){return(x)}                    ##### Define Lebesgue measure   
#'MDEResult <- Koul2StageMde(Y,X, "default", Lx, q, Lx)
#'MDE1stageResult <- MDEResult$MDE1stage
#'MDE2stageResult <- MDEResult$MDE2stage
#'
#'beta1 <- MDE1stageResult$betahat1stage
#'residual1 <- MDE1stageResult$residual1stage
#'rho1 <- MDE1stageResult$rhohat1stage
#'
#'beta2 <- MDE2stageResult$betahat2stage
#'residual2 <- MDE2stageResult$residual2stage
#'rho2 <- MDE2stageResult$rhohat2stage



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

Koul2StageMde <- function(Y,X,D,RegIntMeasure, AR_Order, ArIntMeasure){
  
  
  DimMat <- dim(X)
  n <- DimMat[1]
  p <- DimMat[2]
  
  MDE1Result <- KoulLrMde(Y,X, D, RegIntMeasure)
  beta1 <- MDE1Result$betahat
  resid1 <- MDE1Result$residual
  rho1 <- KoulArMde(resid1, AR_Order, ArIntMeasure)$rhohat
  
  MDE1 <- list(betahat1stage=beta1, residual1stage=resid1, rhohat1stage=rho1)
  
  ###########################   2 stage MDE
  Ytilde <- vector(length=(n-AR_Order))
  Xtilde <- matrix(rep(0,times=(n-AR_Order)*p), nrow=(n-AR_Order), ncol=p )
  
  for(j in 1:(n-AR_Order)){
    
    tempX <- rep(0, times=p)
    tempY <- 0
    for (k in 1: AR_Order){
      tempX <- tempX + rho1[k]*Y[AR_Order+j-k, ]
      tempY <- tempY + rho1[k]*Y[AR_Order+j-k]
    }
    Xtilde[j, ] <- X[(j+AR_Order), ] - tempX
    Ytilde[j] <- Y[j+AR_Order] - tempY
  }
  MDE2Result <- KoulLrMde(Ytilde, Xtilde, D, RegIntMeasure)
  beta2 <- MDE2Result$betahat
  resid2 <- Y-X%*%beta2
  rho2 <- KoulArMde(resid2, AR_Order, ArIntMeasure)$rhohat
  
  MDE2 <- list(betahat2stage=beta2, residual2stage=resid2, rhohat2stage=rho2)
  
  ResultVal <- list(MDE1stage=MDE1, MDE2stage=MDE2)
  return(ResultVal)
  
}




#' Minimum distance estimation in linear regression model.
#'
#' Estimates the regression coefficients in the model \eqn{Y=X\beta + \epsilon}.
#'@param Y - Vector of response variables in linear regression model.
#'@param X - Design matrix of explanatory variables in linear regression model.
#'@param D - Weight matrix. Dimension of D should match that of X. "default" uses XA where A=(X'X)^(-1/2).
#'@param IntMeasure - Symmetric and \eqn{\sigma}-finite measure.
#'@return betahat   - Minimum distance estimator of \eqn{\beta}. 
#'@return residual  - Residuals after minimum distance estimation. 
#'@examples
#'####################
#'n <- 10
#'p <- 3
#'X <- matrix(runif(n*p, 0,50), nrow=n, ncol=p)   #### Generate n-by-p design matrix X 
#'beta <- c(-2, 0.3, 1.5)                        #### Generate true beta = (-2, 0.3, 1.5)' 
#'eps <- rnorm(n, 0,1)                           #### Generate errors from N(0,1)
#'Y <- X%*%beta + eps
#'
#'D <- "default"                                 #### Use the default weight matrix
#'
#'Lx <- function(x){return(x)}                   ##### Define Lebesgue measure 
#'MDEResult <- KoulLrMde(Y,X,D, Lx)              ##### Use Lebesgue measure for the integration
#'betahat <- MDEResult$betahat                   ##### Obtain minimum distance estimator 
#'resid <- MDEResult$residual                    ##### Obtain residual 
#'
#'Dx <- function(x){                             ##### Define degenerate measure at 0
#'        if(x==0){
#'          return(1)
#'        }else{
#'          return(0)
#'        }
#'      }
#'MDEResult <- KoulLrMde(Y,X,D, Dx)              ##### Use degenerate measure at 0 for the integration
#'betahat <- MDEResult$betahat                   ##### Obtain minimum distance estimator 
#'resid <- MDEResult$residual                    ##### Obtain residual

#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@seealso KoulArMde() and Koul2StageMde()
#'@export
#'@importFrom "stats" "optim"
#'@importFrom "stats" "nlm"



KoulLrMde <- function(Y, X, D, IntMeasure){
  
  if(nargs() != 4){
    message("Number of arguments should be four.")
    stop()
  }
  Hx = IntMeasure

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
      for (i in 1:nXRow){
        if(nXCol == 1){xi <- X[i]} else {xi <- t(X[i,])}  
        
        dik <- D[i,k]
        ei <- Y[i] - xi %*% t
        expr_i <- expression(Hx(ei))
        ei <- eval(expr_i)
        for(j in i:nXRow){
          if(nXCol == 1){xj <- X[j]} else {xj <- t(X[j,])}
          
          djk <- D[j,k]
          ej <- Y[j] - xj %*% t
          
          expr_j <- expression(Hx(ej))
          ej <- eval(expr_j)
          
          fval <- fval + 2*dik*djk* ( max(ei, -ej)+max(-ei, ej) - max(ei, ej)-max(-ei, -ej))
          
        }
      }
      
      
    }
    return(fval)
    
  }
  
  return(Dual)
}



#' Minimum distance estimation in the autoregression model of the known order.
#'
#' Estimates the autoressive coefficients in the \eqn{X_t = \rho' Z_t + \xi_t } where \eqn{Z_t} is the vector of \eqn{q} observations at times \eqn{t-1,...,t-q}.  
#'@param X - Vector of \code{n} observed values.
#'@param AR_Order - Order of the autoregression model.
#'@param IntMeasure - Symmetric and \eqn{\sigma}-finite measure.
#'@return rhohat - Minimum distance estimator of \eqn{\rho}.
#'@return residual - Residuals after minimum distance estimation. 
#'@examples
#'##### Generate stationary AR(2) process with 10 observations 
#'n <- 10
#'q <- 2
#'rho <- c(-0.2, 0.8)    ##### Generate true parameters rho = (-0.2, 0.8)'
#'eps <- rnorm(n, 0,1)   ##### Generate innovations from N(0,1)
#'X <- rep(0, times=n)
#'for (i in 1:n){
#'  tempCol <- rep(0, times=q)
#'  for (j in 1:q){
#'    if(i-j<=0){
#'      tempCol[j] <- 0
#'    }else{
#'      tempCol[j] <- X[i-j]
#'    }
#'  } 
#'X[i] <- t(tempCol)%*% rho + eps[i]
#'}
#'
#'Lx <- function(x){return(x)}                   ##### Define Lebesgue measure 
#'MDEResult <- KoulArMde(X, q, Lx)               ##### Use Lebesgue measure for the integration
#'rhohat <- MDEResult$rhohat                     ##### Obtain minimum distance estimator
#'resid  <- MDEResult$residual                   ##### Obtain residual
#'
#'Dx <- function(x){                             ##### Define degenerate measure at 0
#'        if(x==0){
#'          return(1)
#'        }else{
#'          return(0)
#'        }
#'      }
#'MDEResult <- KoulArMde(X, q, Dx)               ##### Use degenerate measure at 0 for the integration
#'rhohat <- MDEResult$rhohat                     ##### Obtain minimum distance estimator 
#'resid <- MDEResult$residual                    ##### Obtain residual


#'@references
#'[1] Koul, H. L (1985). Minimum distance estimation in linear regression with unknown error distributions. Statist. Probab. Lett., 3 1-8.
#'@references
#'[2] Koul, H. L (1986). Minimum distance estimation and goodness-of-fit tests in first-order autoregression. Ann. Statist., 14 1194-1213. 
#'@references
#'[3] Koul, H. L (2002). Weighted empirical process in nonlinear dynamic models. Springer, Berlin, Vol. 166
#'@importFrom "stats" "nlminb"
#'@export
#'@seealso KoulLrMde() and Koul2StageMde()





KoulArMde <- function(X, AR_Order, IntMeasure){
  
  Hx = IntMeasure

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

  rho_hat <- Tmin$par
  resid <- Xres - Xexp%*% rho_hat
  
  lst <- list(rhohat=rho_hat, residual=resid)
  
  return(lst)
}

TARLoss <- function(X, AR_Order, Hx){
  
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
        expr_i <- expression(Hx(ei))
        ei <- eval(expr_i)
        
        
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
          expr_j <- expression(Hx(ej))
          ej <- eval(expr_j)
          
          fval <- fval + 2*dik*djk* ( max(ei, -ej)+max(-ei, ej) - max(ei, ej)-max(-ei, -ej))
        }
        
      }      

    }
    return(fval)
    
  }
  
  return(Dual)
}


