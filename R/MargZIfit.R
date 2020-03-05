#' Fit marginal ZI regression model passed from a formula interface or directly using the model parameters
#' @param bmat Design matrix for the non-degenerate component
#' @param gmat Design matrix for the degenerate component
#' @param y    Response variable
#' @param params Vector of intial parametrs to the model
#' @param method Method for fitting: one of ES_BinMix (Es binary mixing), ES_MultiMix (ES multinomial mixing), RS_MultiMix (RS multinomial mixing)
#' @param mulink Link function for the non-degenerate component ("log","sqrt","logit","probit","identity","cloglog","inverse")
#' @param plink Link function for the degenerate component ("logit","probit")
#' @param formula2 Formula object for the degenerate component
#' @param dims Vector or list of parameter dimensions by component
#' @param id Cluster fic ID variable
#' @param time Time or measurement occassion variable
#' @param nClust Number of clusters (sample size)
#' @param clustSize Number of observations in a cluster
#' @param bin Vector of binomial denominators
#' @param errDist Distribution for the non-degenerate component ("poisson","binomial")
#' @param corstruct Correlation structure for the non-degenerate component ("indep","exch","ar1","unspec","onedep","toep3"), default = "indep"
#' @param mixcorstr Correlation structure for the degenerate component ("indep","exch","ar1","unspec","onedep","toep3"), default = "indep"
#' @param nophi Logical: Whether dispersion parameter is to be estimated or not, default = FALSE
#' @param verbose Logical: Wheteher to print out model fitting steps, default = FALSE
#' @return list of model parameters and standard errors.

margZIfitEst <- function(bmat, gmat, y, params, method, mulink, plink="logit",  formula2,dims, id, time, nClust, clustSize, bin,
                         errDist, corstruct="indep", mixcorstr="indep", nophi, verbose=FALSE){


  # Order of params should be: params <- c(beta, gamma, rho/phi, delta)
  dimbeta <- dims$dimbeta; dimgamma <- dims$dimgamma
  N = nClust*clustSize
  R <- 2^clustSize

  dimrho <- dimdelta <- dimphi <- 0
  if(corstruct!='indep') dimrho <- dims$dimrho
  if(mixcorstr!='indep' & method=="ESBinMix") dimdelta <- dims$dimdelta
  if(!nophi) dimphi <- 1

  #Extract the parameters
  beta <- params[1:dimbeta];  gamma <- params[(dimbeta+1):(dimbeta+dimgamma)]; phi <- 1
  rho <- delta <- 0
  if(corstruct!='indep') rho <- params[(dimbeta+dimgamma+1):(dimbeta+dimgamma+dimrho)]
  if(mixcorstr!='indep' & method=="ESBinMix") delta <- params[(dimbeta+dimgamma+dimrho+dimphi+1):(dimbeta+dimgamma+dimrho+dimphi+dimdelta)]
  if(!nophi) phi <- params[(dimbeta+dimgamma+dimrho+dimphi)]

  dimparm <- dimgamma + dimbeta + dimrho + dimphi + dimdelta
  dimparmB <- dimbeta + dimrho + dimphi;  dimparmG <- dimgamma +dimdelta

  if(method=="ESBinMix"){
    model.result <- zipgeeES1(bmat, gmat, y, params, mulink, plink, dims, id, time, nClust,
                              clustSize, bin, errDist, corstruct, mixcorstr, nophi, verbose=FALSE)
  }else if(method=='ESMultiMix'){

    ## Change gamma parameter and dimensions
    gammainit <- rep(1,R-1); dimgamma <- length(gammainit)
    ## Update the dimensions List
    dims$dimgamma <- dimgamma

    # Initial parameter estimates
    params <- c(beta,gammainit)

    # Update the mode; matrices and weights to dimensions of the data
    bMat <- as.matrix(bmat[1:nClust,]); Bin <- bin[1:nClust];
    gMat <- as.matrix(gmat[1:nClust,1])

    if(corstruct!='indep') params <- c(params,rho)
    if(!nophi) params <- c(params,phi)

    # Compute Initial estimates  for gamma
    Y0 <- as.numeric(y==0)
    Y0Mat <- matrix(Y0, nrow = nClust, ncol = clustSize, byrow = FALSE)
    Y0Tab <- ftable(as.data.frame(Y0Mat), col.vars = 1:clustSize)
    Noclass <- which(Y0Tab==0)
    R <- 2^clustSize
    len_Noclass <- length(Noclass); R2 <- R-len_Noclass

    model.result <- MultiMixES1(bmat, gmat, y, params, mulink, plink, formula2, dims, id, time, nClust,
                                clustSize, bin, errDist, corstruct, nophi, verbose=FALSE)

    }else if(method=="RSMultiMix"){
    ## Update gamma and dimension
    gammainit <- rep(1,R-1); dimgamma <- length(gammainit)

    ## Update the dimensions List
    dims$dimgamma <- dimgamma

    ## Initial parameter estimates
    params <- c(beta,gammainit)
    ## Update the mode; matrices and weights to dimensions of the data
    bMat <- as.matrix(bmat[1:nClust,]); Bin <- bin[1:nClust];
    gMat <- as.matrix(gmat[1:nClust,1])

    if(corstruct!='indep') params <- c(params,rho)
    if(!nophi) params <- c(params,phi)

    ## Compute Initial estimates  for gamma
    BinY0 <- as.numeric(y==0)
    BinY0Mat <- matrix(BinY0, nrow = nClust, ncol = clustSize, byrow = FALSE)
    BinY0Tab <- ftable(as.data.frame(BinY0Mat), col.vars = 1:clustSize)
    Noclass <- which(BinY0Tab==0)
    R <- 2^clustSize
    len_Noclass <- length(Noclass); R2 <- R-len_Noclass

    model.result <- ModelBasedZIfit(bmat, gmat, y, params, mulink, plink, dims, id, time, nClust,
                                     clustSize, bin, errDist, corstruct, nophi, verbose=FALSE)
  }

  ## SUMMARIZE MODEL RESULTS and COMPUTE STANDARD ERRORS
  # number of iterations; Yes/No Convergence
  iterations <- model.result$'iterations'
  convergence <- model.result$'converge'
  beta_fit <- model.result$'Beta'; gamma_fit <- model.result$'Gamma'
  names(beta_fit) <- colnames(bmat)
  if(method=="ESBinMix") names(gamma_fit) <- colnames(gmat)

  # Correlation estimates: Poisson; Degenerate
  corparm_fit <- rho_fit <- delta_fit <- 0;
  phi_fit <- 1
  if(corstruct!='indep'){
    rho_fit <- model.result$'rho'
    names(rho_fit) <- rep("count_correlation",dimrho)
    corparm_fit <- rho_fit
  }
  if(mixcorstr!='indep' & method=="ESBinMix"){
    delta_fit <- model.result$'delta'
    names(delta_fit) <- rep("zero_correlation",dimdelta)
    if(corstruct!='indep'){
      corparm_fit <- c(corparm_fit,delta_fit)
      }else{
        corparm_fit <- delta_fit
      }
  }

  if(!nophi){
    phi_fit <- model.result$'phi'
    names(phi_fit) <- rep("dispersion",dimphi)
  }

  # STANDARD ERRORS
  # Final parameter estimates
  params_fit <- c(beta_fit, gamma_fit)
  if(corstruct!='indep') params_fit <- c(params_fit, rho_fit)
  if(!nophi) params_fit <- c(params_fit, phi_fit)
  if(mixcorstr!='indep' & method=="ESBinMix") params_fit <- c(params_fit, delta_fit)

  if(method=="ESBinMix"){
    # Compute the standard errors
    UList <- list()
    for(i in 1:nClust){
      UList[[i]] <- getScorei(params_fit, index=i, y, dims, mulink, plink, bmat, gmat, id, time, nClust,
                              clustSize, bin, errDist, corstruct, mixcorstr, nophi)
    }
    VmatList <- lapply(UList, function(x) x%*%t(x))
    Vmat <- Reduce("+", VmatList)

    gradList <- list()
    for(i in 1:nClust){
      gradList[[i]] <- numDeriv::jacobian(getScorei, x=params_fit, method="simple", index=i, y=y, dims=dims,
                                mulink=mulink, plink=plink, bmat=bmat, gmat=gmat, id=id,
                                time=time, nClust=nClust, clustSize=clustSize, bin=bin,
                                errDist=errDist, corstruct=corstruct,mixcorstr=mixcorstr, nophi=nophi)
    }
    gradS <- Reduce("+", gradList)
    eta_hat <- solve(gradS)
    VarCovParmsHat <- eta_hat%*%Vmat%*%t(eta_hat)
    colnames(VarCovParmsHat) <- rownames(VarCovParmsHat) <- names(params_fit)
    StdErr_fit <- sqrt(diag(VarCovParmsHat))
    count_StdErr_fit <- StdErr_fit[1:(dimbeta)]
    zero_StdErr_fit <- StdErr_fit[(dimbeta+1):(dimbeta+dimgamma)]
    corparm_StdErr <- 0
    if(corstruct!="indep"){
      rho_StdErr <- StdErr_fit[(dimbeta+dimgamma+1):(dimbeta+dimgamma+dimrho)]
      corparm_StdErr <- rho_StdErr
    }
    if(mixcorstr!="indep"){
      delta_StdErr <- StdErr_fit[(dimbeta+dimrho+dimgamma+1):(dimbeta+dimrho+dimphi+dimgamma+dimdelta)]
      if(corstruct!='indep'){
        corparm_StdErr <- c(corparm_StdErr,delta_StdErr)
        }else{
          corparm_StdErr <- delta_StdErr
        }
    }
  }else if(method=="ESMultiMix"){
    dims$dimgamma <- length(gamma_fit)
    UiList <- list()
    for(i in 1:nClust){
      UiList[[i]] <- getMultiScorei(params_fit, Index=i, y=y, dims=dims,mulink=mulink, plink=plink,
                                    bmat=bmat, gmat=gMat, id=id, time=time, nClust=nClust,
                                    clustSize=clustSize, bin=bin, errDist=errDist,
                                    corstruct=corstruct, nophi=nophi)
    }

    VList <- lapply(UiList, function(x) x%*%t(x))
    VMat <- Reduce("+", VList)

    gradList <- list()
    for(i in 1:nClust){
      gradList[[i]] <- numDeriv::jacobian(getMultiScorei,x = params_fit, method="simple", Index=i, y=y,
                                dims=dims,mulink=mulink, plink=plink, bmat=bmat, gmat=gMat,
                                id=id, time=time, nClust=nClust, clustSize=clustSize,
                                bin=bin, errDist=errDist, corstruct=corstruct, nophi=nophi)
    }
    gradU <- Reduce("+", gradList)
    eta_hat <- solve(gradU)
    VarCovParmsHat <- eta_hat%*%VMat%*%t(eta_hat)
    colnames(VarCovParmsHat) <- rownames(VarCovParmsHat) <- names(params_fit)
    VarParm <- diag(VarCovParmsHat)
    StdErr_fit <- sqrt(VarParm)
    count_StdErr_fit <- StdErr_fit[1:(dimbeta)]
    zero_StdErr_fit <- StdErr_fit[(dimbeta+1):(dimbeta+dims$dimgamma)]
    corparm_StdErr <- 0
    if(corstruct!="indep") corparm_StdErr <- StdErr_fit[(dimbeta+dims$dimgamma+1):(dimbeta+dims$dimgamma+dimrho)]
  }else if(method=="RSMultiMix"){
    dims$dimgamma <- length(gamma_fit)

    omegaU <- getEstep(params_fit, y=y, dims=dims,mulink=mulink, plink=plink,
                       bmat=bmat, gmat=gMat, id=id, time=time, nClust=nClust,
                       clustSize=clustSize, bin=bin, errDist=errDist,
                       corstruct=corstruct, nophi=nophi)

    UiList <- list()
    for(i in 1:nClust){
      UiList[[i]] <- getMultiMBScorei(params_fit, Index=i, omega=omegaU$omega, u=omegaU$u, y=y, dims=dims,
                                      mulink=mulink, plink=plink, bmat=bmat, gmat=gMat, id=id, time=time,
                                      nClust=nClust, clustSize=clustSize, bin=bin, errDist=ERRDIST2,
                                      corstruct=corstruct, nophi=nophi)
    }

    VList <- lapply(UiList, function(x) x%*%t(x))
    VMat <- Reduce("+", VList)

    gradList <- list()
    for(i in 1:nClust){
      gradList[[i]] <- numDeriv::jacobian(getMultiMBScorei,x = params_fit, method="simple", Index=i,omega=omegaU$omega,
                                u=omegaU$u, y=y, dims=dims,mulink=mulink, plink=plink, bmat=bmat, gmat=gMat,
                                id=id, time=time, nClust=nClust, clustSize=clustSize,bin=bin,
                                errDist=errDist, corstruct=corstruct, nophi=nophi)
    }
    gradU <- Reduce("+", gradList)
    eta_hat <- solve(gradU)
    VarCovParmsHat <- eta_hat%*%VMat%*%t(eta_hat)
    colnames(VarCovParmsHat) <- rownames(VarCovParmsHat) <- names(params_fit)
    VarParm <- diag(VarCovParmsHat)
    StdErr_fit <- sqrt(VarParm)
    StdErr_fit <- sqrt(VarParm)
    count_StdErr_fit <- StdErr_fit[1:(dimbeta)]
    zero_StdErr_fit <- StdErr_fit[(dimbeta+1):(dimbeta+dims$dimgamma)]
    corparm_StdErr <- 0
    if(corstruct!="indep") corparm_StdErr <- StdErr_fit[(dimbeta+dims$dimgamma+1):(dimbeta+dims$dimgamma+dimrho)]
  }

list(iterations=iterations,converge=convergence,count_coefficients=beta_fit,
     zero_coefficients=gamma_fit, count_correlation=rho_fit,
     correlation_coefficients=corparm_fit, correlation_StdErr=corparm_StdErr,
     zero_correlation=delta_fit, phi=phi_fit, coefficients=params_fit,
     VarCov=VarCovParmsHat, count_StdErr_fit=count_StdErr_fit,
     zero_StdErr_fit=zero_StdErr_fit)

}

margZIfit <- function(bmat, ...) UseMethod("margZIfit")

margZIfit.default <- function(bmat, gmat, y, params, method, mulink, plink, formula2, dims, id, time, nClust,
                              clustSize, bin, errDist, corstruct, mixcorstr, nophi, verbose=FALSE){


  bmat <- as.matrix(bmat)
  gmat <- as.matrix(gmat)
  y <- as.numeric(y)
  params <- as.numeric(params)
  est <- margZIfitEst(bmat, gmat, y, params, method, mulink, plink,formula2, dims, id, time, nClust, clustSize, bin,
                      errDist, corstruct, mixcorstr, nophi=nophi, verbose=FALSE)
  est$call <- match.call()
  class(est) <- "margZIfit"
  est
}


#' Formula interface for fitting marginal ZI regression models
#'
#' @param formula Two part formula, first part non-degenerate, second part degenerate  e.g. y ~ 1 + a | 1
#' @param data Specify data
#' @param idvar Specify the cluster specific ID variable
#' @param timevar Specify the time or measurement occassion variable
#' @param method Method for fitting: one of ES_BinMix (Es binary mixing), ES_MultiMix (ES multinomial mixing), RS_MultiMix (RS multinomial mixing)
#' @param mulink Link function for the non-degenerate component ("log","sqrt","logit","probit","identity","cloglog","inverse")
#' @param plink Link function for the degenerate component ("logit","probit")
#' @param errDist Distribution for the non-degenerate component ("poisson","binomial")
#' @param corstruct Correlation structure for the non-degenerate component ("indep","exch","ar1","unspec","onedep","toep3"), default = "indep"
#' @param mixcorstr Correlation structure for the degenerate component ("indep","exch","ar1","unspec","onedep","toep3"), default = "indep"
#' @param nophi Logical: Whether dispersion parameter is to be estimated or not, default = FALSE

#'
#' @return model parameters and standard errors
#' @export

margZIfit.formula <- function(formula, data=list(), idvar,timevar, method, mulink, plink="logit", errDist, corstruct="indep",
                              mixcorstr="indep", nophi){

  id <- data[,idvar]
  time <- data[ ,timevar]


  ## Extract the model inputs from formula
  FORMULA <- Formula(formula)
  modframe <- model.frame(FORMULA, data=data)
  ## Extract model matrices
  bmat <- model.matrix(FORMULA, data = modframe, rhs = 1)
  gmat <- model.matrix(FORMULA, data = modframe, rhs = 2)
  y <- model.response(modframe)
  formula1 <- formula(FORMULA, rhs=1)
  ExtractFormula2 <- formula(FORMULA, rhs=2)
  TermsFormula2 <- terms(ExtractFormula2)
  formula2 <- formula(delete.response(TermsFormula2))
  nClust <- length(unique(id)); clustSize <- length(unique(time)); N <- nClust*clustSize; R <- 2^clustSize

  bin <- rep(1,N);

  # Generate Initial Parameters From Independent Models
  if(errDist=="poisson"){
    #IndepModel <- zeroinfl(formula, data = LongData)
    IndepModel <- gamlss::gamlss(formula1, formula2, data = data, family = ZIP(mu.link = "log", sigma.link = "logit"),
                         trace=FALSE)
  }else if(errDist=="binomial"){
    IndepModel <- gamlss::gamlss(formula1, formula2, data = data, family = ZIBI(mu.link = "logit", sigma.link ="logit"),
                         trace=FALSE)
  }

  # The response matrix
  Ymat <- matrix(y, ncol = clustSize , byrow = FALSE)

  # Define dimensions
  dimbeta <- ncol(bmat)
  dimgamma <- ncol(gmat)

  # Set initial parameters for the nondegenerate part
  betainit <- IndepModel$mu.coefficients; rhoinit <- 0.5; phiinit <- 1
  names(betainit) <- colnames(bmat)
  # Set initial parameters for the nondegenerate part
  gammainit <- IndepModel$sigma.coefficients; deltainit <- 0.5
  names(gammainit) <- colnames(gmat)
  # More dimensions
  dimbeta <- length(betainit); dimgamma <- length(gammainit)
  dimdelta <- dimrho <- dimphi <- 0
  if(corstruct!='indep') dimrho <- length(rhoinit)
  if(!nophi) dimphi <- length(phiinit)
  if(mixcorstr!='indep') dimdelta <- length(deltainit)
  dims <- list(dimbeta=dimbeta,dimgamma=dimgamma,dimrho=dimrho,dimdelta=dimdelta,dimphi=dimphi)

  params <- c(betainit, gammainit)
  if(corstruct!='indep') params <- c(params,rhoinit)
  if(!nophi) parm <- c(params,phiinit)
  if(mixcorstr!='indep') params <- c(params,deltainit)

  method <- method

  est <- margZIfit.default(bmat=bmat, gmat=gmat, y=y, params=params, method=method, mulink=mulink, plink=plink,
                           formula2=formula2, dims=dims, id=id, time=time, nClust=nClust, clustSize=clustSize, bin=bin,
                           errDist=errDist, corstruct=corstruct, mixcorstr=mixcorstr, nophi=nophi,verbose=FALSE)
  est$call <- match.call()
  est$formula <- formula
  est
}


summary.margZIfit <- function(object, ...){
  count_params <- object$count_coefficients
  count_StdErr <- object$count_StdErr_fit
  count_zvalue <-  count_params/count_StdErr
  count_pvalue <- 2*pnorm(abs(count_zvalue), lower.tail = FALSE)
  count_CoefTab <- cbind(Estimate = count_params,
                   StdErr = count_StdErr,
                   z.value = count_zvalue,
                   p.value = count_pvalue)
  zero_params <- object$zero_coefficients
  zero_StdErr <- object$zero_StdErr_fit
  zero_zvalue <-  zero_params/zero_StdErr
  zero_pvalue <- 2*pnorm(abs(zero_zvalue), lower.tail = FALSE)
  zero_CoefTab <- cbind(Estimate = zero_params,
                         StdErr = zero_StdErr,
                         z.value = zero_zvalue,
                         p.value = zero_pvalue)
  if(is.null(object$correlation_coefficients)){
    correlation_CoefTab <- NULL
  }else{
    correlation_params <- object$correlation_coefficients
    correlation_StdErr <- object$correlation_StdErr
    correlation_zvalue <-  correlation_params/correlation_StdErr
    correlation_pvalue <- 2*pnorm(abs(correlation_zvalue), lower.tail = FALSE)
    correlation_CoefTab <- cbind(Estimate = correlation_params,
                                 StdErr = correlation_StdErr)
    #z.value = correlation_zvalue,
    #p.value = correlation_pvalue)
    }
  dispersion <- object$phi
  convergence <- object$converge
  iterations <- object$iterations
  result <- list(call=object$call, count_coefficients=count_CoefTab,
                 zero_coefficients=zero_CoefTab,
  correlation_coefficients=correlation_CoefTab, dispersion=dispersion,
  convergence=convergence, iterations=iterations)
  class(result) <- "summary.margZIfit"
  result
}

print.summary.margZIfit <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Count model coefficients\n")
  cat("\n")
  printCoefmat(x$count_coefficients, P.value=TRUE, has.Pvalue=TRUE)
  cat("\n")
  cat("Zero-Inflation model coefficients \n")
  cat("\n")
  printCoefmat(x$zero_coefficients, P.value=TRUE, has.Pvalue=TRUE)
  cat("\n")
  if(!is.null(x$correlation_coefficients)){
    cat("Correlation coefficients \n")
    cat("\n")
    printCoefmat(x$correlation_coefficients, P.value=FALSE, has.Pvalue=FALSE)
    }#else{
  #   x$correlation_coefficients <- NULL
  # }
  cat("\n")
  cat("The number of iterations:",x$iterations)
  cat("\n")
  cat("Convergence:", x$convergence)
  cat("\n")
  cat("Dispersion parameter taken to be",x$dispersion)
  cat("\n")
}

