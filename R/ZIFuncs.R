# ZIFuncs.R

# This script contains functions that are used in several different algorithms to fit
# marginal ZI regression models of various forms that arise in Tawanda's dissertation.

# Better to have one file where all these functions reside than have the functions in different
# scripts. This way there should just be one version of each of these functions which avoids
# having to maintain multiple versions or worry about differences between versions.

# The functions here are:

# getmu(x, bin=1, link = "no default") - calculates mean vector of a GLM with given link,
#                                        binomial denominator bin, and linear predictor x

# getmixp(x,link=PLINK) - calculates mixing prob vector from linear predictor x under a given
#                         link

# getf2(y,mu,bin=1,errDist=ERRDIST2) - calculates density of y for a given mean mu and binomial
#                                      denominator for given error distn. Used for getting
#                                      density of y under the non-degenerate distribution in the
#                                      mixture model.

# getvmu(mu,errDist=ERRDIST2) - calculates variance function for a given mean mu under a given
#                               error distribution

# getdvmu(mu,errDist=ERRDIST2) - calculates derivative of variance function

# getdmu(x,link="no default") - calculates derivative of mean for a GLM with respect to the
#                               linear predictor, evaluates at linear predictor x

# getri(alpha,ni,corstruct=corstruct) - calculates working correlation matrix for a cluster
#                                       of size ni evaluated at a correlation parameter alpha

# getdvidalpha(ni,alpha,vmui,alphaindex,phi,corstruct=corstruct)
#   - Calculates derivative of cluster-specific variance matrix (vi) with respect to alpha

# GausCor(Sig) - I believe this function computes V_i22 in the notation of Hall's (2001)
#                paper on EQL from V_i11 which is passed as the matrix Sig to this function.
#                This GausCor() function seems not to be used at all elsewhere in our functions.
#                It is not a necessary function to implement EGEE, because EGEE can be written
#                as equation (15) in Hall (2001), which does not involve V_i22

# vecnorm(x) - computes sqrt(sum(x^2)), the norm of the vector x

# expit(x) - computes the inverse logit or expit of x

# vech(mat) - computes the vech of a square matrix mat.
#             It stackes the columns of the lower triangle of mat.

# vec(mat) - computes the vec of a square matrix mat. It stacks the columns of mat.

# getScoreHessi(dta, corParm, dispParm = 1,nophi = nophi,errDist = ERRDIST2,link = MU2LINK,
#               corstruct = corstruct,weighted=FALSE)
#   - This function computes the contribution to the EGEE estimating function for a single
#     cluster, which I am calling efi, and the contribution to the expected gradient of the
#     score vector (plays role of Hessian) for the ith cluster, which I am calling ii.
#     This function will be lapply'ed to a list of data frames
#     Each data frame in this list should have the following variables:
#         ID = subject/cluster ID
#         timeIndex = time index giving ordering of observations within each cluster
#         X1-Xp = cols of design matrix (in order), including an intercept as X1, if desired.
#         Y = variable playing the role of the response
#         EY = variable playing role of the expected response
#         lp = linear predictor
#         u = conditional mixing weights

# posSeq(a,b) - a simple little function to compute a sequence from a to b, returning a vector
#               of length 0 if b<a. Needed in mvpoisrng().

# mvpoisrng(sigma) - Generates a random poisson vector with var-cov matrix sigma using the
#                    method of Sim (reference?) This function works but has not been revised to
#                    improve the code (e.g., for efficiency)

# getSigma(lam, rho) - Currently unused, but had been used to implement Sim's method for
#                      correlated Poissons

# revlogit(linpred) - I don't know what this function is for and it seems unused

# getmultip(x) - This function is included here but commented out. It seems to be
#                incorrectly written.

# getMultip(lpmat) - A corrected version of getmultip(). Computes generalized expit to give
#             matrix of category probabilities from a matrix of linear predictors corresponding
#             to a generalized logit model for an R category response. lpmat has R-1 columns
#             and so does matrix that this function returns





getmu = function(x, bin=1, link = "no default"){
  #Calculates mean vector of a GLM with link link given a linear predictor x
  # bin is the binomial denom
  switch(link,logit=bin/(1+exp(-x)),probit=bin*pnorm(x),log=bin*exp(x),identity=bin*x,
         inverse=bin/x,cloglog=bin*(1-exp(-exp(x))),
         "invalid link in getmu()")
}

getmixp = function(x,link){
  #computes mixing probability from linear predictor based on the chosen link
  #function for the mixing probability (PLINK)
  #x = GMAT%*%gam

  switch(link,logit=1/(1+exp(-x)),probit=pnorm(x),"invalid link in getmixp()")
}

getf2 = function(y,mu,bin=1,errDist){

  switch(errDist,poisson=dpois(y,mu),binomial=dbinom(y,bin,mu/bin),
         "Invalid error distribution in getf2()")
}

getvmu = function(mu,errDist){
  #Calculates variance function from mean mu.

  switch(errDist,normal=rep(1,length(mu)),poisson=mu,binomial=mu*(1-mu),
         gamma=mu*mu,"invalid error distribution in getvmu()")
}

getdmu = function(x,link="no default"){
  #Calculates derivative of mean vector of a GLM with link link with respect to linear
  # predictor x

  switch(link,logit=exp(x)/(1+exp(x))^2,log=exp(x),identity=rep(1,length(x)),inverse=-1/x/x,
         "invalid link in getdmu()")
}

getdvmu = function(mu,errDist){
  # Calculates derivative of variance function.
  # where mu is the mean vector

  switch(errDist,normal=rep(0,length(mu)),poisson=rep(1,length(mu)),binomial=1-2*mu,
         gamma=2*mu,"invalid error distribution in getdvmu()")
}

getri = function(alpha,ni,corstruct){
  #Calculates correlation matrix of current GLM for level 2 unit i
  #USE:getri(alpha,ni)
  #where alpha is the association vector, ni is ith cluster size

  m <- diag(ni)
  if(corstruct == 'unspec'){
    ritran <-  m
    ritran[lower.tri(ritran)]=alpha
    ri <- ritran + t(ritran) - m
  } else if(corstruct == 'toep3'){
    dimalpha = length(alpha)
    if(dimalpha != 3 | ni != 4)  stop('dim of alpha must be 3 and cluster size must equal 4 to use toep3')
    ri <- matrix( c(1,alpha)[c(1+abs(row(m)-col(m)))],nrow=4)
  } else {
    ri <- switch(corstruct,exch=alpha*(1-m)+m,ar1=alpha^abs(row(m)-col(m)),
                 onedep=m+alpha*((abs(row(m)-col(m))==1)+0),indep=m,
                 "invalid corstruct in getri()")
  }
  return(ri)
}

getdvidalpha = function(ni,alpha,vmui,alphaindex,phi,corstruct=corstruct){
  # Calculates derivative of variance matrix with respect to alpha

  dimalpha = length(alpha)
  m <- diag(ni)

  if(alphaindex<=dimalpha){
    if(corstruct=='exch'){
      dri = matrix(1,ni,ni)-diag(ni)
    } else if(corstruct=='indep'){
      dri = diag(ni)
    } else  if(corstruct=='ar1'){
      dri <- abs(row(m)-col(m))* alpha^(abs(row(m)-col(m))-1)
    } else if(corstruct=='onedep'){
      dri <- ((abs(row(m)-col(m))==1)+0)
    } else if(corstruct=='unspec'){
      dri = matrix(0,ni,ni)
      indtran <-  dri
      indtran[lower.tri(indtran)]=1:length(alpha)
      indmat <- indtran + t(indtran)
      dri[indmat==alphaindex] <- 1
    }  else if(corstruct=='toep3'){
      dri = matrix(0,ni,ni)
      dri[abs(row(dri)-col(dri))==alphaindex] <- 1
    }

    aisqrt=diag(sqrt(as.vector(vmui)))
    dvi=phi*aisqrt%*%dri%*%aisqrt
  } else if(alphaindex>dimalpha){

    ri=getri(alpha,ni,corstruct)
    aisqrt=diag(sqrt(as.vector(vmui)))
    dvi=aisqrt%*%ri%*%aisqrt
  }
  return(dvi)
}

#Gaussian Correlation Structure - we seem not to use this function as far as I can tell
GausCor = function(Sig){ # Sig=V_i11 passed as a matrix
  ni = nrow(Sig)
  q=ni*(ni+1)/2
  M=matrix(-1, nrow = ni,ncol = ni)
  cov = matrix(0,nrow=q,ncol=q)
  M[lower.tri(M,diag=TRUE)] <- 1:q
  # M= matrix of indicators (need to code this)
  for(x in 1:q){
    j=which(M==x,arr.ind=T)[1]
    k=which(M==x,arr.ind=T)[2]
    for(y in 1:q){
      l=which(M==y,arr.ind=T)[1]
      m=which(M==y,arr.ind=T)[2]

      cov[x,y] = Sig[j,l]*Sig[k,m]+Sig[j,m]*Sig[k,l]
    }
  }
  return(cov)
}

vecnorm = function(x) sqrt(sum(x^2))

expit <- function(x){1/(1+exp(-x))}

vech <- function(mat){
  #returns vech of mat. mat must be square and is assumed symmetric
  # stacks columns of lower triangle of mat

  if(nrow(mat) != ncol(mat)) stop("mat not square in vech()")
  mat[lower.tri(mat,diag=TRUE)]
}

vec = function(mat){
  #returns vec of mat. mat must be square and is assumed symmetric

  #if(nrow(mat) != ncol(mat)) stop("mat not square in vec()")
  mat[1:dim(mat)[1]^2]
}

# Need a function that computes the contribution to the estimating
# function for a single cluster, which I am calling efi, and the
# contribution to the expected gradient of the score vector (plays role of
# Hessian) for the ith cluster, which I am calling ii.
# Function will be lapply'ed to a list of data frames
# Each data frame in this list should have the following variables:
# ID = subject/cluster ID
# timeIndex = time index giving ordering of observations within each cluster
# X1-Xp = columns of the design matrix (in order), including an intercept as X1, if desired.
# Y = variable playing the role of the response
# EY = variable playing role of the expected response
# lp = linear predictor
# u = conditional mixing weights
getScoreHessi <- function(dta, corParm, dispParm = 1, nophi, errDist, link,
                          corstruct, weighted=FALSE){
  #weighti <- dta[,grepl("weight",colnames(dta))][1,]
  ni <- length(dta$id) # cluster size
  q <- ni*ni # dimension of quadratic part of the estimating function
  #q <- ni*(ni+1)/2
  vech=vec # note that here we use a formulation of EGEE involving vec() not vech()

  mui <- dta$EY                         # mean response for ith cluster had been called pi
  lpi <- dta$lp                         # linear predictor for ith cluster
  if (weighted) ui <- dta$u            # Conditional mixing weights
  Yi <- dta$y                           # formerly ui
  Xi <- as.matrix(dta[,substr(names(dta),1,1)=="X"])  #Covariates

  dimMeanParm <- ncol(as.matrix(Xi))
  dimCorParm <- length(corParm)
  dimDispParm <- as.numeric(!nophi)
  dimParm <- dimMeanParm+dimCorParm+dimDispParm

  vmui <- getvmu(mui,errDist)           #had been called vpi
  ri <- getri(corParm,ni,corstruct)
  aisqrt <- diag(sqrt(vmui))
  vi <- dispParm*aisqrt%*%ri%*%aisqrt
  viinv <- solve(vi)
  f1i <- Yi-mui
  ti <- vech(f1i%*%t(f1i));
  taui <- vech(vi)
  f2i <- ti-taui
  fi <- c(f1i, f2i)
  wi <- rbind(cbind(viinv, matrix(0,nrow = ni, ncol = q)),
              cbind(matrix(0, nrow = q, ncol = ni), diag(q)))
  dmui <- getdmu(lpi,link)
  dvmui <- getdvmu(mui,errDist)
  dmuidparmt <- diag(as.vector(dmui))%*%Xi # ni x dim of regression parm
  di22 <- matrix(0, nrow = q, ncol = (dimParm-dimMeanParm))
  ti22 <- di22

  #compute dvi w.r.t. corParm(j)
  for(j in 1:dimCorParm){
    dvij <- getdvidalpha(ni,corParm,vmui,j,dispParm,corstruct)
    dviinvj <- -viinv%*%dvij%*%viinv
    di22[,j] <- vech(dviinvj)
    ti22[,j] <- vech(dvij)
  }
  if(!nophi){
    di22[,(dimCorParm+dimDispParm)] <- vech(-(1/dispParm)*viinv)
    ti22[,(dimCorParm+dimDispParm)] <- vech((1/dispParm)*vi)
  }

  ti12 <- matrix(0, nrow = q, ncol = dimMeanParm)
  for(k in 1:dimMeanParm){
    dvidparm <- dispParm*(diag(0.5*vmui^(-0.5)*dvmui*dmuidparmt[,k])%*%ri%*%aisqrt +
                            aisqrt%*%ri%*%diag(0.5*vmui^(-0.5)*dvmui*dmuidparmt[,k]))
    ti12[,k] <- vech(dvidparm)
  }
  di <- rbind(cbind(dmuidparmt, matrix(0, nrow = ni, ncol = (dimParm-dimMeanParm))),
              cbind(matrix(0, nrow = q, ncol = dimMeanParm), di22))

  ei <- rbind(cbind(dmuidparmt, matrix(0,nrow = ni, ncol = (dimParm-dimMeanParm))),
              cbind(ti12, ti22))


  if(!weighted){
    #weights1 <- vec(vi+Yi%*%t(Yi)-Yi%*%t(mui)-mui%*%t(Yi)+mui%*%t(mui))
    weights1 <- vech((Yi-mui)%*%t(Yi-mui))
    delta_part <- ExpS_4 <- Term_4 <- 0
    if(corstruct!='indep'){
      gamma_part <- t(dmuidparmt)%*%viinv%*%f1i
      ExpS_2 <- t(dmuidparmt)%*%viinv%*%mui
      delta_part <- t(di22)%*%f2i
      ExpS_4 <- t(di22)%*%weights1
      Term_4 <- t(di22)%*%taui
      Term_2 <- t(dmuidparmt)%*%viinv%*%Yi
    } else {
      viinv2 <- solve(dispParm*aisqrt%*%aisqrt)
      gamma_part <- t(dmuidparmt)%*%viinv2%*%f1i
      ExpS_2 <- t(dmuidparmt)%*%viinv2%*%mui
      Term_2 <- t(dmuidparmt)%*%viinv2%*%Yi
    }

    efi <- t(di)%*%wi%*%fi
    ii <- -t(di)%*%wi%*%ei
  } else {
    hi1 <- diag(1-ui)
    #if(MB) {
      #hi2 <- diag(weighti) } else {
    hi2 <- diag(vech((1-ui)%*%t(1-ui)))
    #}
    hi <- diag( c(diag(hi1),diag(hi2))  )
    #weights2 <- vec((1-ui)%*%t(1-ui))
    rho_part <- ExpS_3 <- Term_3 <- 0
    if(corstruct!='indep'){
      ExpS_1 <- -t(dmuidparmt)%*%viinv%*%diag(ui)%*%f1i
      beta_part <- t(dmuidparmt)%*%viinv%*%hi1%*%f1i
      ExpS_3 <- t(di22)%*%hi2%*%f2i
      rho_part <- ExpS_3
      Term_1 <- -t(dmuidparmt)%*%viinv%*%f1i
      Term_3 <- rep(0, dimCorParm)
    } else {
      viinv2 <- solve(dispParm*aisqrt%*%aisqrt)
      ExpS_1 <- -t(dmuidparmt)%*%viinv2%*%diag(ui)%*%f1i
      Term_1 <- -t(dmuidparmt)%*%viinv2%*%f1i
      beta_part <- t(dmuidparmt)%*%viinv2%*%hi1%*%f1i
    }

    efi <- t(di)%*%wi%*%hi%*%fi
    ii <- -t(di)%*%wi%*%hi%*%ei
  }

  if(!weighted){
   out <- list("efi"=efi,"ii"=ii,"ExpS_2"=ExpS_2,"ExpS_4"=ExpS_4,"Term_2"=Term_2,"Term_4"=Term_4,
               "gamma_part"=gamma_part,"delta_part"=delta_part)
  }else{

   out <- list("efi"=efi,"ii"=ii,"ExpS_1"=ExpS_1,"ExpS_3"=ExpS_3,"Term_1"=Term_1,"Term_3"=Term_3,
               "beta_part"=beta_part,"rho_part"=rho_part)
  }
  return(out)
}

wtSum <- function(omega=omegaFull, indMat=jkMat, b2Mat=b2Mat){
  wtvec <- numeric(nrow(indMat))
  for (ind in 1:nrow(indMat)){
    j <- indMat[ind,1]; k <- indMat[ind,2]
    wtvec[ind] <- sum(omega[ b2Mat[,j]==0 & b2Mat[,k]==0])
  }
  wtvec
}

posSeq <- function(a,b){seq(from=a,to=b,length.out=max(0,b+1-a))}

mvpoisrng <- function(sigma){
  #% Generates a random poisson vector with var-cov matrix sigma using the method
  #% of Sim

  if(dim(sigma)[1]!=dim(sigma)[2]) {
    stop('Invalid var-cov matrix in mvpoisrng.m')

  }
  p = dim(sigma)[1]

  lambda = rep(0, p)

  x = rep(0, p)
  z = rep(0, p)

  alpha = matrix(0, nrow = p, ncol = p)

  lambda[1] = sigma[1, 1]

  for (j in 2:p) {
    alpha[j, 1] = sigma[1, j] / lambda[1]

    if (alpha[j, 1] < 0 || alpha[j, 1] > 1) {
      stop('alpha_j1 out of bounds in mvpoisrng.m')

    }
    for (i in posSeq(2,j-1)) {
      a = 0

      for (k in posSeq(1,i-1)) {
        a = a + alpha[i, k] * alpha[j, k] * lambda[k]

      }
      cat("alpha is\n")
      print(alpha)
      cat("a is\n")
      print(a)
      if ((sigma[i, j] > a) && (sigma[i, j] < (lambda[i] + a))) {
        alpha[j, i] = (sigma[i, j] - a) / lambda[i]
      } else{
        stop('alpha_ji not in [0,1] in mvpoisrng.m')

      }
    }
    b = 0

    for (k in posSeq(1,j-1)) {
      b = b + alpha[j, k] * lambda[k]

    }
    if (sigma[j, j] > b) {
      lambda[j] = sigma[j, j] - b
    } else{
      stop('lambda_i is non-positive in mvpoisrng.m')

    }
  } #%j
  x[1] = rpois(1, lambda[1])

  z[1]=x[1]

  for (j in 1:p) {
    sum = 0

    for (i in posSeq(1,j - 1)) {
      if (x[i] > 0) {
        sum = sum + rbinom(1, x[i], alpha[j, i])

      }
    }
    x[j] = rpois(1, lambda[j])

    z[j] = sum + x[j]

  }

  ranpoisvec = z
  return(ranpoisvec)
}

getSigma = function(lam, rho){
  # Needs a comment explaining this function's purpose

  A = diag(length(lam))
  A[lower.tri(A)] = sqrt(prod(lam))*rho
  A[upper.tri(A)] = A[lower.tri(A)]
  return(A)
}

revlogit = function(linpred){
  # Needs a comment explaining this function's purpose

  glin = apply(linpred,2,expit)
  l = apply(glin,1,function(x) 1-sum(x))
  glin = cbind(glin,l)
  return(glin)
}

getmultip = function(x) exp(x)/(1+sum(exp(x)))
# This function looks wrong, should it be like so?:

getMultip <- function(lpmat) {
  # computes inverse generalized logit fom a matrix of linear predictors lpmat
  # Columns of lpmat are assumed to be the R-1 linear predictors for a generalized logit model
  # for response with R categories. Rows are assumed to correspond to different observations
  exp(lpmat)/(1+rowSums(exp(lpmat)))
}

getMultiGMAT <- function(gmat, nClust, index){
  l = ncol(gmat)
  Xi = matrix(0, nrow = R2-1, ncol = l*(R2-1))
  count = 1
  for(j in 1:(R2-1)){
    Xi[j,count:(j*l)] <-gmat[index,]
    count = 1 + j*l
  }
  return(Xi)
}

getMultiGMAT2 <- function(gmat, nClust, index){
  l = ncol(gmat)
  Xi = matrix(0, nrow = R2, ncol = l*R2)
  count = 1
  for(j in 1:R2){
    Xi[j,count:(j*l)] <-gmat[index,]
    count = 1 + j*l
  }
  return(Xi)
}
#ID: Subject ID
#glp: vectorized linear predictor
#y: vectorized response
#EYvec: vectorized probabbilites
MultinomFuncs <- function(dta,weighted=TRUE){
  glpVec <- dta$glpVec
  y <- dta$y
  EYVec <- dta$EYvec
  Xi <- as.matrix(dta[,substr(names(dta),1,1)=="X"])  #Covariates
  vmu <- EYVec*(1-EYVec)
  Vmat <- -outer(EYVec,EYVec)
  diag(Vmat) <- vmu
  invVmat <- solve(Vmat)
  Ddbeta <- Vmat
  weight <- Ddbeta%*%invVmat%*%t(Ddbeta)
  score = Xi%*%Ddbeta%*%invVmat%*%(y-EYVec)
  termDat <- Xi%*%Ddbeta%*%invVmat%*%y
  termParm <- Xi%*%Ddbeta%*%invVmat%*%EYVec
  fish = score%*%t(score)
  dscore = -Xi%*%weight%*%t(Xi)
  return(list(efdot=score, fisher = fish, idot=dscore,termParm=termParm,termDat=termDat))
}

# SigCode <- function(pv){
#   if(pv >= 0 & pv <= 0.001){
#     sig <- '***'}
#   else if(pv > 0.001 & pv <= 0.01){
#     sig <- '**'}
#   else if(pv > 0.01 & pv <= 0.05){
#     sig <- '*'}
#   else if(pv > 0.05 & pv <= 0.1){
#     sig <- '.'}
#   else{ sig <- ' '}
#   return(sig)
# }


# SumParm <- function(params, StdErr){
#
#   cat('Call:' , StrFormula1,'\n')
#   zVal <-  params/StdErr
#   pVal <- 2*pnorm(abs(zVal), lower.tail = FALSE)
#   cat('Coefficients:\n')
#   TabSum <- round(data.frame(cbind(params, StdErr, zVal, pVal)),4)
#   TabSum <- cbind(TabSum,unlist(lapply(pVal,SigCode)))
#   colnames(TabSum) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)','')
#   print(TabSum1 <- TabSum[1:dimbeta,])
#   cat('Signif. codes:', "0 '***'" , "0.001 '**'", "0.01 '*'", "0.05 '.'", "0.1 ' '", "1\n")
#   cat('------\n')
#   cat('Call:' , StrFormula2,'\n')
#   print(TabSum1 <- TabSum[(dimbeta+1):nrow(TabSum),])
#   cat('Signif. codes:', "0 '***'" , "0.001 '**'", "0.01 '*'", "0.05 '.'", "0.1 ' '", "1\n")
# }
