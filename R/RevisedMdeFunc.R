#' Minimum distance estimation in linear regression model.
#'
#' Estimates the regression coefficients in the model \eqn{Y=X\beta + \epsilon}.
#'@param Y - Vector of response variables in linear regression model.
#'@param X - Design matrix of explanatory variables in linear regression model.
#'@param D - Weight Matrix. Dimension of D should match that of X. Default value is XA where A=(X'X)^(-1/2).
#'@param b0 - Initial value for beta.
#'@param IntMeasure - Symmetric and \eqn{\sigma}-finite measure.
#'@return betahat   - Minimum distance estimator of \eqn{\beta}.
#'@return residual  - Residuals after minimum distance estimation.
#'@examples
#'####################
#'n <- 10
#'p <- 3
#'X <- matrix(runif(n*p, 0,50), nrow=n, ncol=p)  #### Generate n-by-p design matrix X
#'beta <- c(-2, 0.3, 1.5)                        #### Generate true beta = (-2, 0.3, 1.5)'
#'eps <- rnorm(n, 0,1)                           #### Generate errors from N(0,1)
#'Y <- X%*%beta + eps
#'
#'D <- "default"                                 #### Use the default weight matrix
#'b0 <- solve(t(X)%*%X)%*%(t(X)%*%Y)             #### Set initial value for beta
#'Lx <- "Lebesgue"                               ##### Define Lebesgue measure
#'MDEResult <- KoulLrMde(Y,X,D, b0, Lx)          ##### Use Lebesgue measure for the integration
#'betahat <- MDEResult$betahat                   ##### Obtain minimum distance estimator
#'resid <- MDEResult$residual                    ##### Obtain residual
#'
#'Dx <- "Degenerate"                             ##### Define degenerate measure at 0
#'
#'MDEResult <- KoulLrMde(Y,X,D, b0, Dx)          ##### Use degenerate measure for the integration
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
#'@importFrom Rcpp evalCpp
#'@useDynLib KoulMde



KoulLrMde <- function(Y, X, D, b0, IntMeasure){

  if(nargs() != 5){
    message("Number of arguments should be five.")
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

    if (nXCol != length(b0) ){
      message("b0 is not conformable to X.")
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

  iter <- 1000
  critVal <- 0.001

  if(IntMeasure == "Lebesgue"){
    bhat <- EstimateBetaMDESimple(Y, X, D, b0, iter, critVal)
  }else if(IntMeasure == "Degenerate"){
    bhat <- EstimateDegenBetaMDE(Y, X, D, b0, iter, critVal)
  }else{
    message("Integrating measure should be either Lebesgue or Degenerate.")
    stop()
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





FindPlace <- function(E, al){

  OriginalE <- E
  LenE <- length(E)

  if(al < E[1]){

    E <- c(E[1]-1, E[1])
    IndexSNum <- 0
    IndexENum <- 1

    lst <- list(E, IndexSNum, IndexENum)
    return(lst)

  }

  if(al > E[LenE]){

    E <- c(E[LenE], E[LenE]+1)

    IndexSNum <- LenE
    IndexENum <- LenE+1

    lst <- list(E, IndexSNum, IndexENum)
    return(lst)

  }




  IndexSNum <- 1
  IndexENum <- LenE

  trialnum <- floor(log(LenE, 2) )+1


  for(i in 1:trialnum){
    LenE <- floor(LenE/2)

    if(al < E[LenE]){

      #lstL <- GetLowerInterval(E,LenE, IndexSNum, IndexENum)

      lstL <- GLI(E,LenE, IndexSNum, IndexENum)

      E <- lstL[[1]]
      IndexSNum <- lstL[[2]]
      IndexENum <- lstL[[3]]
      LenE <- IndexENum-IndexSNum

    }else{
      #lst <- GetUpperInterval(E,LenE, IndexSNum, IndexENum)
      lst <- GUI(E,LenE, IndexSNum, IndexENum)
      E <- lst[[1]]
      IndexSNum <- lst[[2]]
      IndexENum <- lst[[3]]
      LenE <- IndexENum-IndexSNum


      if(E[1]>al){
        IndexENum <- IndexSNum
        IndexSNum <- IndexSNum-1

        E <- OriginalE[IndexSNum:IndexENum]
        break
      }

    }

    if( (IndexENum-IndexSNum) == 1){break}
  }
  lst <- list(E, IndexSNum, IndexENum)
  return(lst)


}






FindIndex <- function(IndLen, EPlus, EMinus, PointVec, gVec, hVec){

  SIndex <- 1
  EIndex <- IndLen

  SVal <- 0
  EVal <- 0


  OldSlopeVal <- 0

  trialnum <- ceiling(log(IndLen, 2) )+1

  #print(trialnum)

  for(i in 1:trialnum){

    MiddleIndex <- ceiling((SIndex+EIndex)/2)

    glst <- FindPlace(EPlus, PointVec[MiddleIndex])

    #glst <- FP(EPlus, PointVec[MiddleIndex])

    #print(glst[[1]])
    #print(glst[[2]])
    #print(glst[[3]])

    #stop()


    if(EMinus == 0){
      SlopeVal <- gVec[( glst[[3]]) ]
    }else{
      hlst <- FindPlace(EMinus, PointVec[MiddleIndex])
      #hlst <- FP(EMinus, PointVec[MiddleIndex])
      SlopeVal <- gVec[( glst[[3]]) ] - hVec[(hlst[[3]])]
    }

    if(SlopeVal > 0){
      EIndex <- MiddleIndex
      EVal <- SlopeVal

    }else{
      SIndex <- MiddleIndex
      Sval <- SlopeVal
    }

    if( (EIndex-SIndex) == 1){break}
  }

  return(SIndex)

}


GetbVec <- function(Y,X, DYM, XpmM, bVec, l){

  dimX = dim(X)
  n=dimX[1]
  p=dimX[2]


  Dstar = DYM[(1:n), ]
  LossVal = CLoss(Y,X,Dstar, bVec)


  pmm = PMM(n, DYM, XpmM, (l-1), p, bVec)

  PMLen = pmm[(2*n^2+1),1]
  MMLen = pmm[(2*n^2+1),2]


  PMat = pmm[ ( 1:PMLen), ]
  if(MMLen != 0){
    MMat = pmm[ ((n^2+1):(n^2+MMLen)), ]
  }

  E = rep(0, times=(PMLen+MMLen))
  E[1:PMLen] = PMat[1:PMLen,3]
  if(MMLen != 0){E[ (PMLen+1) : (PMLen+MMLen) ] = MMat[1:MMLen,3]}


  SortedE = sort(E)
  LenE = length(E)



  PointVec = PV(SortedE, LenE)

  PMat = PMat[order(PMat[,3]),]
  EPlus = PMat[,3]
  if(MMLen != 0){
    MMat = MMat[order(MMat[,3]),]
    EMinus = MMat[,3]
  }else{
    EMinus=0
  }


  gVec = GVec(PMat)

  if(MMLen != 0){
    hVec = HVec(MMat)
  }





  SlopeVec = rep(0,times = (LenE+1))

  ########
  ValVec = rep(0,times = (LenE+1))
  ##############


  TriedbVec = bVec
  nstar = FindIndex((LenE+1), EPlus, EMinus, PointVec, gVec, hVec)

  if( nstar > 0){
    TriedbVec[l] = SortedE[nstar]
    TriedLossVal = CLoss(Y,X,Dstar,TriedbVec)

    Lmat2Val = 0

    if(TriedLossVal < LossVal ){

      Lmat2Val = TriedLossVal-LossVal
      LossVal = TriedLossVal
      bVec=TriedbVec
    }
  }

  lst=list()
  lst[[1]] = bVec
  lst[[2]] = Lmat2Val
  return(lst)


}

EstimateBetaMDESimple <- function(Y, X, D, b0, iter, critVal){
  dimX = dim(X)
  n = dimX[1]
  p = dimX[2]


  DYM = DY(Y, D)
  XpmM = Xpm(X)

  bVec = b0

  LMat = matrix(rep(0, times=p*2), nrow=p, ncol=2)
  for(i in 1:p){LMat[i,1]=i}

  for(ith in 1:iter){
    OldbVec = bVec

    if(ith == 1){

      for(lprime in 1 : p){
        l = LMat[lprime,1]

        ResultLst =  GetbVec(Y,X, DYM, XpmM, bVec, l)
        bVec = ResultLst[[1]]
        LMat[l,2] = ResultLst[[2]]
      }
    }else{
      l = LMat[1,1]

      ResultLst =  GetbVec(Y,X, DYM, XpmM, bVec, l)
      bVec = ResultLst[[1]]
      LMat[1,2] = ResultLst[[2]]
    }


    LMat=  LMat[order(LMat[,2]),]

    dv = bVec-OldbVec
    normdiff = sum(abs(dv)^2)^(1/2)

    #normdiff = norm( bVec-OldbVec , type="2")

    if(normdiff < critVal){break}

  }

  return(bVec)

}



#################################
GetDegenbVec <- function(Y,X, D, bVec, l){

  dimX = dim(X)
  n=dimX[1]
  p=dimX[2]

  LossVal = TLRLossDegen(Y,X,D)(bVec)
  nElLen = 0
  ElMat = matrix(rep(0, times=(2*n)), nrow=n, ncol=2 )

  for(i in 1:n){
    xil = X[i,l]
    if(xil != 0){
      nElLen = nElLen+1
      ElMat[nElLen,1]=nElLen
      ElMat[nElLen,2]= (Y[i]- X[i,]%*%bVec + xil*bVec[l])/xil
    }
  }
  TriedbVec = bVec

  TriedLossVal = LossVal
  for(i in 1:nElLen){
    TriedbVec[l] = ElMat[i,2]
    TestVal = TLRLossDegen(Y,X,D)(TriedbVec)

    if( TestVal < TriedLossVal ){
      bVec[l] = ElMat[i,2]
      TriedLossVal = TestVal
    }

  }

  return(bVec)
}



EstimateDegenBetaMDE <- function(Y, X, D, b0, iter, critVal){

  dimX = dim(X)
  n = dimX[1]
  p = dimX[2]

  tempA = (t(X)%*%X)
  A = sqrtmat(tempA, -0.5)
  D = X%*%A


  bVec = b0

  for(ith in 1:iter){
    OldbVec = bVec
    for(l in 1:p){
      bVec =  GetDegenbVec(Y,X, D, bVec, l)
    }
    dv = bVec-OldbVec
    #normdiff = norm(c, type="2")
    normdiff = sum(abs(dv)^2)^(1/2)

    if(normdiff < critVal){break}

  }


  return(bVec)

}


TLRLossDegen <- function(Y, X, D){

  if (is.vector(X) == TRUE ){
    nXRow <- length(X)
    nXCol <- 1

  }else {
    DimMat <- dim(X)
    LengY <- length(Y)

    nXRow <- DimMat[1]
    nXCol <- DimMat[2]

  }

  p <- nXCol

  Dual <- function(t){
    fval <- 0

    for(k in 1:p){
      tempVal <- 0
      for (i in 1:nXRow){
        if(nXCol == 1){xi <- X[i]} else {xi <- t(X[i,])}
        ei <- Y[i] - xi %*% t
        if(ei > 0){
          ei <- 1
        }else if(ei == 0){
          ei <- 0
        }else{
          ei <- (-1)
        }
        tempVal <- tempVal + D[i,k]*ei
      }
      fval <- fval + tempVal^2
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
#'Lx <- "Lebesgue"                               ##### Define Lebesgue measure
#'MDEResult <- KoulArMde(X, q, Lx)               ##### Use Lebesgue measure for the integration
#'rhohat <- MDEResult$rhohat                     ##### Obtain minimum distance estimator
#'resid  <- MDEResult$residual                   ##### Obtain residual
#'
#'Dx <- "Degenerate"                             ##### Define degenerate measure at 0
#'MDEResult <- KoulArMde(X, q, Dx)               ##### Use degenerate measure for the integration
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

  if ( (Hx != "Lebesgue") && (Hx != "Degenerate") ){
    message("Integrating measure should be either Lebesgue or Degenerate.")
    stop()
  }

  nLength <- length(X)

  if(nLength <= AR_Order){
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
        if(Hx == "Degenerate"){ei <- abs(ei)}


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
          if(Hx == "Degenerate"){ej <- abs(ej)}

          if(j == i){
            fac=1
          }else{
            fac=2
          }

          fval <- fval + fac*dik*djk* ( max(ei, -ej)+max(-ei, ej) - max(ei, ej)-max(-ei, -ej))
        }

      }

    }

    return(fval)

  }

  return(Dual)
}



#'Two-stage minimum distance estimation in linear regression model with autoregressive error.
#'
#'Estimates both regression and autoregressive coefficients in the model \eqn{Y=X\beta + \epsilon} where \eqn{\epsilon} is autoregressive process of known order \code{q}
#'@param Y - Vector of response variables in linear regression model.
#'@param X - Design matrix of explanatory variables in linear regression model.
#'@param D - Weight Matrix. Dimension of D should match that of X. Default value is XA where A=(X'X)^(-1/2).
#'@param b0 - Initial value for beta.
#'@param RegIntMeasure - Symmetric and \eqn{\sigma}-finite measure used for estimating \eqn{\beta}. Either Lebesgue or Degenerate measure.
#'@param AR_Order - Order of the autoregressive error.
#'@param ArIntMeasure - Symmetric and \eqn{\sigma}-finite measure used for estimating autoregressive coefficients of the error. Either Lebesgue or Degenerate measure.
#'@return MDE1stage - The list of the first stage minimum distance estimation result. It contains betahat1stage, residual1stage, and rho1stage.
#'\itemize{
#'  \item betahat1stage - The first stage minimum distance estimators of regression coefficients.
#'  \item residual1stage - Residuals after the first stage minimum distance estimation.
#'  \item rho1stage - The first stage minimum distance estimators of autoregressive coefficients of the error.
#'}
#'@return MDE2stage - The list of the second stage minimum distance estimation result. It contains betahat2stage, residual2stage, and rho2stage.
#'\itemize{
#'  \item betahat2stage - The second stage minimum distance estimators of regression coefficients.
#'  \item residual2stage - Residuals after the second stage minimum distance estimation.
#'  \item rho2stage - The second stage minimum distance estimators of autoregressive coefficients of the error.
#'}
#'@examples
#'####################
#'n <- 10
#'p <- 3
#'X <- matrix(runif(n*p, 0,50), nrow=n, ncol=p)  #### Generate n-by-p design matrix X
#'beta <- c(-2, 0.3, 1.5)                        #### Generate true beta = (-2, 0.3, 1.5)'
#'rho  <- 0.4                                    #### True rho = 0.4
#'eps <- vector(length=n)
#'xi <- rnorm(n, 0,1)                            #### Generate innovation from N(0,1)
#'                                               #### Generate autoregressive process of order 1
#'for(i in 1:n){
#'  if(i==1){eps[i] <- xi[i]}
#'  else{eps[i] <- rho*eps[i-1] + xi[i]}
#'}
#'Y <- X%*%beta + eps
#'#####################
#'D <- "default"                                  #### Use the default weight matrix
#'b0 <- solve(t(X)%*%X)%*%(t(X)%*%Y)              #### Set initial value for beta
#'
#'Lx <- "Lebesgue"                                ##### Define Lebesgue measure
#'MDEResult <- Koul2StageMde(Y,X, "default", b0, Lx, 1, Lx)
#'MDE1stageResult <- MDEResult[[1]]
#'MDE2stageResult <- MDEResult[[2]]
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

Koul2StageMde <- function(Y,X,D, b0, RegIntMeasure, AR_Order, ArIntMeasure){


  DimMat <- dim(X)
  n <- DimMat[1]
  p <- DimMat[2]

  MDE1Result <- KoulLrMde(Y,X, D, b0, RegIntMeasure)
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
      tempX <- tempX + rho1[k]*X[AR_Order+j-k, ]
      tempY <- tempY + rho1[k]*Y[AR_Order+j-k]
    }
    Xtilde[j, ] <- X[(j+AR_Order), ] - tempX
    Ytilde[j] <- Y[j+AR_Order] - tempY
  }
  MDE2Result <- KoulLrMde(Ytilde, Xtilde, D, beta1, RegIntMeasure)
  beta2 <- MDE2Result$betahat
  resid2 <- Y-X%*%beta2
  rho2 <- KoulArMde(resid2, AR_Order, ArIntMeasure)$rhohat

  MDE2 <- list(betahat2stage=beta2, residual2stage=resid2, rhohat2stage=rho2)

  ResultVal <- list(MDE1stage=MDE1, MDE2stage=MDE2)
  return(ResultVal)

}









