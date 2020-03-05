
zipgeeES1 <- function(bmat, gmat, y, params, mulink, plink, dims, id, time, nClust,
                      clustSize, bin, errDist, corstruct, mixcorstr, nophi, verbose=FALSE){

  #Control for model convergence
  BIG <- 1.0e10; SMALL <- 1.0e-10
  maxits <- 500; hardmax <- 500; toler <- 1.0e-6;
  converge <- FALSE

  #Order of params should be: params <- c(beta, gamma, rho/phi, delta)
  dimbeta <- dims$dimbeta; dimgamma <- dims$dimgamma
  N = nClust*clustSize

  dimrho <- dimphi <- dimdelta <- 0
  if(corstruct!='indep') dimrho <- dims$dimrho
  if(mixcorstr!='indep') dimdelta <- dims$dimdelta
  if(!nophi) dimphi <- 1

  #Extract the parameters
  beta <- params[1:dimbeta];  gamma <- params[(dimbeta+1):(dimbeta+dimgamma)]; phi <- 1
  rho <- delta <- 0
  if(corstruct!='indep') rho <- params[(dimbeta+dimgamma+1):(dimbeta+dimgamma+dimrho)]
  if(!nophi) phi <- params[(dimbeta+dimgamma+dimrho+dimphi)]
  if(mixcorstr!='indep'&method=="ESBinMix") delta <- params[(dimbeta+dimgamma+dimrho+dimphi+1):(dimbeta+dimgamma+dimrho+dimphi+dimdelta)]

  dimparm <- dimgamma + dimbeta + dimrho + dimphi + dimdelta
  dimparmB <- dimbeta + dimrho + dimphi;  dimparmG <- dimgamma + dimdelta

  parm <- gamma
  if(mixcorstr!='indep') parm <- c(parm,delta)
  parm <- c(parm,beta)
  if(corstruct!='indep') parm <- c(parm,rho)
  if(!nophi) parm <- c(parm,phi)

  # Now start the loop that iterates the E and S steps
  its <- 0
  while(!converge){
    its <- its+1
    if (verbose) cat("ES iteration number ",its,"\n")
    if (verbose) (parmold <- parm)
    parmold <- parm
    #Compute E(Z|Y) The E-step weiights
    blp <- bmat%*%beta; glp <- gmat%*%gamma
    p <- getmixp(glp, link = plink)
    mu2 <- getmu(blp, bin = bin, link = mulink)
    u <- rep(0,N)
    fofyij <- 1 + ((1-p)/p)*getf2(y, mu = mu2, bin = bin, errDist = errDist)
    u[y==0] <- 1/fofyij[y==0]


    #M-step for gamma and rho
    if(mixcorstr=='indep'){
      glmfit <- suppressWarnings(glm.fit(gmat, u, family = binomial(link = plink)))
      gamma <- glmfit$coef

    } else {

      geeits <- 0
      geeConv <- FALSE
      while (geeits < maxits) {
        geeits <- geeits + 1
        glp <- gmat %*% gamma
        p <- getmixp(glp, link=plink)
        dta <- data.frame(id = id, time = time, y = u, EY = p, X = gmat, lp = glp)
        dtaList <- split(dta, dta$id)
        # set up the list of data frames and arguments to getScoreHessi
        iList <- lapply(dtaList, getScoreHessi, corParm = delta, dispParm = 1, nophi = TRUE, errDist = "binomial",
                        link = plink, corstruct = mixcorstr)
        i2List <- unlist(iList, recursive = FALSE)
        odd1 <- grep(".efi", names(i2List))
        efdotList <- i2List[odd1]
        efdot <- Reduce("+", efdotList)
        odd2 <- grep(".ii", names(i2List))
        idotList <- i2List[odd2]
        idot <- Reduce("+", idotList)
        idotinv <- solve(idot)
        updatevec <- idotinv %*% efdot

        maxscore <- max(abs(efdot))
        maxchange <- max(abs(updatevec))

        gamma <- gamma - updatevec[1:dimgamma]
        delta <- delta - updatevec[(dimgamma + 1):(dimgamma + dimdelta)]

        if (verbose) cat("Gamma step ",geeits,"gamma=",gamma,", delta=",delta,"\n")

        maxparm <- max(c(abs(gamma), abs(delta)))
        minparm <- min(c(abs(gamma), abs(delta)))

        if (maxparm > BIG ||minparm < SMALL || rcond(idot) < SMALL) {
          #not converging
          geeits <- maxits + 1
        }
        if (maxscore < toler && maxchange < toler) {
          geeConv <- TRUE
          geeits <- maxits + 1
        }

      } #End while loop for gamma,delta
      if (!geeConv) {
        stop("GEE iterations failed in S step for gamma")
      } else
        geeConv <- FALSE
    } #End else group after if(mixcorstr=='indep')

    #M step for beta, rho, phi
    if(corstruct == 'indep') {
      glmwts = 1 - u
      if(errDist=='poisson'){
        famfun <- if(nophi) poisson(link=mulink) else quasipoisson(link=mulink)
      } else if(errDist=='binomial'){
        famfun <- if(nophi) binomial(link=mulink) else quasibinomial(link=mulink)
      }
      #famfun <- switch(errDist, poisson = ifelse(nophi, poisson(link=mulink), quasipoisson(link=mulink)),
      #                           binomial = ifelse(nophi, binomial(link=mulink), quasibinomial(link=mulink)))
      glmfit <- suppressWarnings(glm.fit(bmat, y, family = famfun,
                                         weights = glmwts))
      beta <- glmfit$coef
      if(!nophi)  phi = summary(glmfit)$dispersion

    } else {
      geeits <- 0
      while (geeits < maxits) {
        geeits <- geeits + 1
        blp <- bmat%*%beta
        mu <- getmu(blp, bin = bin, link = mulink)
        dta2 <- data.frame(id = id, time = time, y = y , EY = mu, X = bmat, lp = blp, u = u)
        dtaList2 <- split(dta2, dta2$id)
        # set up the list of data frames and arguments to getScoreHessi
        iList.2 <- lapply(dtaList2, getScoreHessi, corParm = rho, dispParm = phi, nophi = nophi,
                          errDist = errDist, link = mulink, corstruct = corstruct, weighted = TRUE)

        i2List.2 <- unlist(iList.2, recursive = FALSE)
        odd1 <- grep(".efi", names(i2List.2))
        udotList <- i2List.2[odd1]
        udot <- Reduce("+", udotList)
        odd2 <- grep(".ii", names(i2List.2))
        idotList <- i2List.2[odd2]
        idot <- Reduce("+", idotList)
        idotinv <- solve(idot)
        updatevec <- idotinv %*% udot

        maxscore <- max(abs(udot))
        maxchange <- max(abs(updatevec))

        beta <- beta - updatevec[1:dimbeta]
        rho <- rho - updatevec[(dimbeta + 1):(dimbeta + dimrho)]
        phi <- phi - updatevec[dimparmB] * as.numeric(!nophi)

        if(verbose) cat("Beta step ",geeits,"beta=",beta,", rho=",rho,"phi"=phi,"\n")

        maxparm <- max(c(abs(beta), abs(rho)))
        minparm <- min(c(abs(beta), abs(rho)))

        if (maxparm > BIG || minparm < SMALL || rcond(idot) < SMALL) {
          # not converging
          geeits <- maxits + 1
        }

        if (maxscore < toler && maxchange < toler) {
          geeConv <- TRUE
          geeits <- maxits + 1
        }  #if(!geefail)
      } #End while loop for beta,rho,phi

      if (!geeConv) {
        converge <- FALSE
        stop("S step for beta did not converge")
        #      its = maxits+1
      }
    } #End else group after if(corstruct=='indep')

    parm <- gamma

    if(mixcorstr!='indep') parm <- c(parm,delta)
    parm <- c(parm,beta)
    if(corstruct != 'indep') parm <- c(parm,rho)
    if(!nophi) parm <- c(parm,phi)


    maxchange <- max(abs(parm-parmold))
    relchange <- vecnorm(parm-parmold)/vecnorm(parm)

    if(relchange < toler){
      converge <- TRUE
    }
  } #End global while loop
  if(verbose) cat("Gamma estimates"=gamma,"Beta estimates"=beta,"rho"=rho,"phi"=phi,"delta"=delta,
                  "iterations"=its,"converge"=converge,"\n")

  return(list("Gamma"=gamma,"Beta"=beta, "rho" = rho, "phi" = phi, "delta" = delta,
              "iterations"=its,"converge"=converge))
} #End function


MultiMixES1 <- function(bmat, gmat, y, params, mulink, plink, formula2, dims, id, time, nClust,
                        clustSize, bin, errDist, corstruct, nophi=TRUE, verbose=FALSE){

  BIG <- 1.0e10; SMALL <- 1.0e-10
  #Control for model convergence
  maxits <- 500; hardmax <- 500; toler <- 1.0e-6;
  converge <- FALSE

  #Order of params should be: params <- c(beta, gamma, rho/phi, delta)
  dimbeta <- dims$dimbeta; dimgamma <- dims$dimgamma
  N = nClust*clustSize
  R <- 2^clustSize

  dimrho <- dimphi <- 0
  if(corstruct!='indep') dimrho <- dims$dimrho
  if(!nophi) dimphi <- 1

  #Extract the parameters
  beta <- params[1:dimbeta];  gamma <- params[(dimbeta+1):(dimbeta+dimgamma)]; phi <- 1
  rho <- 0
  if(corstruct!='indep') rho <- params[(dimbeta+dimgamma+1):(dimbeta+dimgamma+dimrho)]
  if(!nophi) phi <- params[(dimbeta+dimgamma+dimrho+dimphi)]

  dimparm <- dimgamma + dimbeta + dimrho + dimphi
  dimparmB <- dimbeta + dimrho + dimphi;  dimparmG <- dimgamma

  # set up initial values for the nu_ir's. This is done by counting how many times each pattern
  # of zero response values happens in the within-cluster data over all clusters. Thus
  # the rows of nuInit are identical (assume initially that nu_ir=nu_r)
  y0 <- as.numeric(y==0)
  y0Mat <- matrix(y0, nrow = nClust, ncol = clustSize, byrow = FALSE)
  y0Tab <- ftable(as.data.frame(y0Mat), col.vars = 1:clustSize)
  Noclass <- which(y0Tab==0)
  len_Noclass <- length(Noclass)
  if(len_Noclass==0) Noclass <- 1000*R
  R2 <- R-len_Noclass
  y0Tab <- y0Tab[,-Noclass]
  nuFullInit <- matrix(c(y0Tab/sum(y0Tab)),nrow = nClust, ncol = R2, byrow = TRUE)
  colnames(nuFullInit) <- R.utils::intToBin(0:(R-1))[-Noclass] # 1 denotes an excess zero
  nuInit <- nuFullInit[,1:(R2-1)]
  gamma <- gamma[-Noclass]

  # set up initial values for the omega_ir's. These values are indicators for whether the
  # rth pattern of zero valued responses has occurred in the ith cluster.
  # so omega_ir is not constant over i
  yi0Tab <- ftable(data.frame(y0Mat, sub = 1:nClust),col.vars = 1:clustSize, row.vars = "sub")[,-Noclass]
  omegaFullInit <- as.matrix(yi0Tab)
  omegaInit <- omegaFullInit[,1:(R2-1)]

  #GMatList <- lapply(1:nClust, function(x) getMultiGMAT(gmat,nClust,x))
  #GMat <- Reduce(rbind,GMatList, accumulate = FALSE)
  GMatList <- lapply(1:nClust, function(x) kronecker(gmat[x,],diag(R2-1)))
  GMat <- Reduce(rbind,GMatList, accumulate = FALSE)


  #Initial weights to the multinomial model
  #Compute the linear predictor for the non-degenerate distribution
  # blp = bmat%*%beta
  # #Compute the mean for the non-degenerate component
  # muInit <- getmu(blp, bin = bin, link = mulink)
  # f2ofy <- getf2(y, muInit, bin = bin, errDist = errDist)
  # f2ofyMat <- matrix(f2ofy, nrow = nClust, ncol = clustSize, byrow = FALSE)
  #
  # EstepTermsInit <- matrix(0, nrow = nClust, ncol = R2)
  # cellIDs <- intToBin(0:(R-1))[-Noclass]
  # for(j in 1:R2){
  #   indices <- unlist(gregexpr("0",cellIDs[j]))
  #   if(all(indices==-1)) indices <- -(1:clustSize)
  #   EstepTermsInit[,j] <- nuFullInit[,j,drop=FALSE]*rowProds(f2ofyMat[,indices,drop=FALSE])*
  #     rowProds(y0Mat[,-indices,drop=FALSE])
  # }
  # omegaFullStart <- EstepTermsInit/rowSums(EstepTermsInit)
  # wts_init <- matrixcalc::vec(omegaFullStart)
  # nu2_init <- kronecker(diag(R2),rep(1,nClust))
  # dfraMe = cbind.data.frame(z = apply(nu2_init, 1, function(x) which(x==1)), wts_init)
  # initData <- mlogit::mlogit.data(dfraMe, shape="wide", sep="",choice="z",alt.levels=1:R2)
  # initfit <- mlogit(z ~ 0|formula2RHS, weights = wts, data = initData )
  # gamma <- initfit$coefficients
  # dims$dimgamma <- length(gamma)

  parm <- c(gamma,beta)
  if(corstruct!='indep') parm <- c(parm,rho)
  if(!nophi) parm <- c(parm,phi)

  # Now start the loop that iterates the E and S steps
  its <- 0

  #GMatList <- lapply(1:nClust, function(x) getMultiGMAT(gmat,nClust,x))
  #GMat <- Reduce(rbind,GMatList, accumulate = FALSE)
  GMatList <- lapply(1:nClust, function(x) kronecker(diag(R2-1),gmat[x,]))
  GMat <- Reduce(rbind,GMatList, accumulate = FALSE)

  while(!converge){

    its <- its+1
    if (verbose) cat("ES iteration number ",its,"\n")
    parmold <- parm

    ##Compute the E-step for multinomial mixing
    #Compute the linear predictor for the non-degenerate distribution
    blp = bmat%*%beta
    #Compute the mean for the non-degenerate component
    mu <- getmu(blp, bin = bin, link = mulink)

    #Estimates of the conditional mixing probs
    if(its==1) nu = nuFullInit #set nu to initial estimates

    #initialize the estimates of the conditional mixing probs
    #omega_est = matrix(0, nrow = K, ncol = R)

    #Compute the conditional mixing probs (E-Step:)
    f2ofy <- getf2(y, mu, bin = bin, errDist = errDist)
    f2ofyMat <- matrix(f2ofy, nrow = nClust, ncol = clustSize, byrow = FALSE)

    EstepTerms <- matrix(0, nrow = nClust, ncol = R2)
    cellIDs <- R.utils::intToBin(0:(R-1))[-Noclass]
    for(j in 1:R2){
      indices <- unlist(gregexpr("0",cellIDs[j]))
      if(all(indices==-1)) indices <- -(1:clustSize)
      EstepTerms[,j] <- nu[,j,drop=FALSE]*rowProds(f2ofyMat[,indices,drop=FALSE])*
        rowProds(y0Mat[,-indices,drop=FALSE])
    }
    omegaFull <- EstepTerms/rowSums(EstepTerms)

    uMat <- matrix(0, nrow = nClust, ncol = clustSize)
    for(t in 1:clustSize) {
      uMat[,t] <- rowSums(omegaFull[,substr(cellIDs,t,t)=="1"])
    }

    u <- matrixcalc::vec(uMat)
    wts <- as.vector(omegaFull)
    nu2 <- kronecker(diag(R2),rep(1,nClust))
    dfM = cbind.data.frame(z = apply(nu2, 1, function(x) which(x==1)), wts)
    DM <- mlogit::mlogit.data(dfM, shape="wide", sep="",choice="z",alt.levels=1:R2)
    fit3 <- mlogit::mlogit(z ~ -1|1, weights = wts, data = DM)
    gamma <- fit3$coefficients
    glpVec <- GMat%*%gamma
    glpMat <- matrix(glpVec, nrow = nClust, ncol = R2-1, byrow = TRUE)

    #Compute the mean for mixing mechanism
    piMat <- getMultip(glpMat)
    piMatFull <- cbind(piMat,1-rowSums(piMat))
    colnames(piMatFull) <- cellIDs # 1 denotes an excess zero
    nu <- piMatFull


    #M step for beta, rho, phi
    if(corstruct == 'indep') {
      glmwts = 1 - u
      if(errDist=='poisson'){
        famfun <- if(nophi) poisson(link=mulink) else quasipoisson(link=mulink)
      } else if(errDist=='binomial'){
        famfun <- if(nophi)  binomial(link=mulink) else quasibinomial(link=mulink)
      }
      glmfit <- suppressWarnings(glm.fit(bmat, y , family = famfun, weights = glmwts))
      beta <- glmfit$coef
      if(!nophi) phi = summary(glmfit)$dispersion
    } else {
      geeits <- 0
      geeConv <- FALSE

      while (geeits < maxits) {
        geeits <- geeits + 1
        blp <- bmat %*% beta
        mu <- getmu(blp, bin = bin, link = mulink)

        dta2 <- data.frame(id = id, time = time, y = y, EY = mu,
                           X = bmat, lp = blp, u = u)
        dtaList2 <- split(dta2, dta2$id)
        # set up the list of data frames and arguments to getScoreHessi
        iList.2 <- lapply(dtaList2, getScoreHessi, corParm = rho, dispParm = phi,
                          nophi = nophi, errDist = errDist, link = mulink, corstruct = corstruct,
                          weighted = TRUE)

        i2List.2 <- unlist(iList.2, recursive = FALSE)
        udotInd <- grep(".efi", names(i2List.2))
        udotList <- i2List.2[udotInd]
        udot <- Reduce("+", udotList)
        idotInd <- grep(".ii", names(i2List.2))
        idotList <- i2List.2[idotInd]
        idot <- Reduce("+", idotList)

        idotinv <- solve(idot)
        updatevec <- idotinv %*% udot

        maxscore <- max(abs(udot))
        maxchange <- max(abs(updatevec))

        beta <- beta - updatevec[1:dimbeta]
        rho <- rho - updatevec[(dimbeta + 1):(dimbeta + dimrho)]
        phi <- phi - updatevec[dimparmB] * as.numeric(!nophi)
        if (verbose) cat("Beta step ",geeits,"beta=",beta,", rho=",rho,"\n")

        maxparm <- max(c(abs(beta), abs(rho)))
        minparm <- min(c(abs(beta), abs(rho)))

        if (maxparm > BIG || minparm < SMALL || rcond(idot) < SMALL) {
          # not converging
          geeits <- maxits + 1
        }

        if (maxscore < toler && maxchange < toler) {
          geeConv <- TRUE
          geeits <- maxits + 1
        }  #if(!geefail)
      } #End while loop for beta,rho,phi

      if(!geeConv) {
        converge <- FALSE
        stop("S step for beta did not converge")
        #      its = maxits+1
      }
    }
    parm = c(gamma,beta)
    if(corstruct != 'indep'){
      parm <- c(parm,rho)
    }

    if(!nophi){
      parm <- c(parm,phi)
    }

    maxchange <- max(abs(parm-parmold))
    relchange <- vecnorm(parm-parmold)/vecnorm(parm)

    if(relchange < toler){
      converge <- TRUE
    }

} #End global while loop

  return(list("Gamma"=gamma,"Beta"=beta,"rho"=rho,"phi"=phi,
              "iterations"=its,"converge"=converge))
} #End function


getRow <- function(indx,Dat,longDat,nGrps=R2){
  rowDat <- Dat[Dat$subj==indx,]
  supDat <- longDat[longDat$subj==indx,]
  n_row <- nrow(rowDat)
  new_Dat <- data.frame(do.call(rbind,lapply(split(rowDat,1:n_row),
                                             function(x) apply(x, 2, function(y) rep(y, n_row)))))
  new_Dat$choice <- as.vector(t(diag(n_row)))
  new_Dat$scenario <- rep(1:n_row, n_row)
  #if(n_row>1) new_Dat$wgts <- ifelse(new_Dat$choice==1, 0.9900,0.0100)
  new_Dat$wgts <- rep(rowDat$wgts, n_row)
  new_Dat <- new_Dat[order(new_Dat$subj,new_Dat$scenario),]
  return(new_Dat)
}

getdr <- function(r,y,t){
  lmat <- gtools::permutations(n=2,r=t,v=c(0,1),repeats.allowed=TRUE)
  lvec <- lmat[r,]
  if (r==1) dr=1
  else{
    dr <- prod(y[as.logical(lvec)]==0)
  }
  dr
}

ModelBasedZIfit <- function(bmat, gmat, y, params, mulink, plink, dims, id, time, nClust,
                           clustSize, bin, errDist, corstruct, nophi=TRUE, verbose=FALSE){

  # Control for model convergence
  BIG <- 1.0e10; SMALL <- 1.0e-10
  maxits = 100;  hardmax = 500; toler = 1.0e-6;
  converge <- FALSE
  #Order of params should be: params <- c(beta, gamma, rho/phi)
  dimbeta <- dims$dimbeta; dimgamma <- dims$dimgamma
  N = nClust*clustSize
  R <- 2^clustSize

  dimrho <- dimphi <- 0
  if(corstruct!='indep') dimrho <- dims$dimrho
  if(!nophi) dimphi <- 1

  #Extract the parameters
  beta <- params[1:dimbeta];  gamma <- params[(dimbeta+1):(dimbeta+dimgamma)]; phi <- 1
  rho <- 0
  if(corstruct!='indep') rho <- params[(dimbeta+dimgamma+1):(dimbeta+dimgamma+dimrho)]
  if(!nophi) phi <- params[(dimbeta+dimgamma+dimrho+dimphi)]

  dimparm <- dimgamma + dimbeta + dimrho + dimphi
  dimparmB <- dimbeta + dimrho + dimphi;  dimparmG <- dimgamma

  # set up initial values for the nu_ir's. This is done by counting how many times each pattern
  # of zero response values happens in the within-cluster data over all clusters. Thus
  # the rows of nuInit are identical (assume initially that nu_ir=nu_r)
  y0 <- as.numeric(y==0)
  y0Mat <- matrix(y0, nrow = nClust, ncol = clustSize, byrow = FALSE)
  y0Tab <- ftable(as.data.frame(y0Mat), col.vars = 1:clustSize)
  Noclass <- which(y0Tab==0)
  len_Noclass <- length(Noclass)
  if(len_Noclass==0) Noclass <- R*1000
  dimgamma <- dimgamma - len_Noclass
  gamma <- gamma[-Noclass]
  cellIDs <- R.utils::intToBin(0:(R-1))[-Noclass] # 1 denotes an excess zero
  R2 <- R-len_Noclass
  y0Tab <- y0Tab[-Noclass]
  nuFullInit <- matrix(c(y0Tab/sum(y0Tab)),nrow = nClust, ncol = R2, byrow = TRUE)
  colnames(nuFullInit) <- cellIDs
  nuInit <- nuFullInit[,1:(R2-1)]


  # set up initial values for the omega_ir's. These values are indicators for whether the
  # rth pattern of zero valued responses has occurred in the ith cluster.
  # so omega_ir is not constant over i
  yi0Tab <- ftable(data.frame(y0Mat, sub = 1:nClust),col.vars = 1:clustSize, row.vars = "sub")[,-Noclass]
  omegaFullInit <- as.matrix(yi0Tab)
  omegaInit <- omegaFullInit[,1:(R2-1)]

  parm <- c(gamma,beta)
  if(corstruct!='indep') parm <- c(parm,rho)
  if(!nophi) parm <- c(parm,phi)

  GMatList <- lapply(1:nClust, function(x) kronecker(diag(R2-1),gmat[x,]))
  GMat <- Reduce(rbind,GMatList, accumulate = FALSE)
  # The response matrix
  Ymat <- matrix(y, ncol = clustSize , byrow = FALSE)
  # Now start the loop that iterates the E and S steps
  its <- 0
  while(!converge){
    its <- its+1
    if(verbose) cat("ES iteration number ",its,"\n")
    parmold <- parm

    # Compute the linear predictor for the non-degenerate distribution
    blp <- bmat%*%beta
    # Compute the mean for the non-degenerate component
    mu <- getmu(blp, bin = bin, link = mulink)

    # Estimates of the conditional mixing probs
    if(its==1){
      nu <- nuFullInit #set nu to initial estimates
      offst <- log(nuFullInit[,1:(ncol(nuFullInit))]/nuFullInit[,ncol(nuFullInit)])
      omega = omegaFullInit
    }

    # Setup the means to compute the covariates
    muMat <- matrix(mu, nrow = nClust , ncol = clustSize)
    #if(verbose) dim(muMat)
    Mulist <- split(muMat, 1:nrow(muMat))
    all_pos <- lapply(strsplit(cellIDs,''),function(a) 1-as.numeric(a))
    all_posList <- rep(all_pos, times=nClust)
    Vechall_pos <- lapply(all_posList, function(x) vech(x%*%t(x)))
    Ylist <- split(Ymat, 1:nrow(Ymat))
    riList <- purrr::map2(Ylist, Mulist, function(x,y) x-y)
    riList_rep <- rep(riList,each=R2)

    # CrossProduct Terms

    # cross_riList <- lapply(riList, function(x) combn(x,m=2, FUN=prod))

    # cross_riList_rep <- rep(cross_riList,each=R2)
    #  #riList2 <- purrr::map2(riList_rep,cross_riList_rep,c)
    suppressWarnings(X_irList <- purrr::map2(riList_rep, Vechall_pos, function(x,y) x*y))

    X_ir <- do.call(rbind,X_irList)


    # Quadratic terms

    VechSqdiff_Yi_Mui <- purrr::map2(Ylist, Mulist, function(x,y) vech((x-y)%*%t(x-y)))
    VechSqdiff_Yi_Mui_Rep <- rep(VechSqdiff_Yi_Mui,each=R2)
    suppressWarnings(X_irListSq <- purrr::map2(VechSqdiff_Yi_Mui_Rep, Vechall_pos, function(x,y) x*y))
    X_irSq <- do.call(rbind,X_irListSq)

    XMat <- cbind(X_ir, X_irSq)

    choice <- as.vector(t(yi0Tab))
    altern <- kronecker(rep(1,nClust),1:R2)
    offst <- as.vector(t(offst))
    subj <- rep(unique(id), each=R2)
    alt_choiceList <- lapply(Ylist, function(x) unlist(lapply(1:R2, getdr, y=x, t=clustSize)))
    alt_choiceMat <- do.call(rbind,alt_choiceList)
    alt_choice <- as.vector(t(do.call(rbind,alt_choiceList)))

    if(its==1){
      #Compute the conditional mixing probs (E-Step:)
      f2ofy <- getf2(y, mu, bin = bin, errDist = errDist)
      f2ofyMat <- matrix(f2ofy, nrow = nClust, ncol = clustSize, byrow = FALSE)

      EstepTerms <- matrix(0, nrow = nClust, ncol = R2)
      cellIDs <- R.utils::intToBin(0:(R-1))[-Noclass]
      for(j in 1:R2){
        indices <- unlist(gregexpr("0",cellIDs[j]))
        if(all(indices==-1)) indices <- -(1:clustSize)
        EstepTerms[,j] <- nu[,j,drop=FALSE]*rowProds(f2ofyMat[,indices,drop=FALSE])*
          rowProds(y0Mat[,-indices,drop=FALSE])
      }
      omega <- EstepTerms/rowSums(EstepTerms)
      wgts <- matrixcalc::vec(t(omega))
    }

    colnames(offst) <- NULL
    ml_dat_long <- data.frame(subj=subj, altern=altern, alt_ch=alt_choice, wgts=wgts,
                              offst=offst,choice=choice, X=XMat)
    row.names(ml_dat_long) <- NULL
    ml_dat_trim <- ml_dat_long[ml_dat_long$alt_ch>0,]
    nalts <- ddply(ml_dat_trim,.variables=~subj,.fun=summarise,
                   nalts=length(subj)) # number of alternatives for each subject
    ml_dat_trim <- merge(ml_dat_trim,nalts,by="subj",all=TRUE)
    row.names(ml_dat_trim) <- NULL

    ml_dat1 <- do.call(rbind,lapply(1:nClust,getRow, Dat=ml_dat_trim,longDat=ml_dat_long,nGrps=R2))
    row.names(ml_dat1) <- NULL
    newsubj <- as.character(ml_dat1[,"subj"])
    newsubjlen <- nchar(newsubj)
    ml_dat1$clust <- paste(newsubj,ml_dat1$scenario ,sep = ".")
    DM_2 <- mlogit::mlogit.data(ml_dat1, shape="long",chid.var="clust",choice = "choice",alt.var="altern")

    PredictorVariables <- paste(c("1",colnames(DM_2)[grepl('X',colnames(DM_2))],"offst"),
                                sep="")
    Formula <- formula(paste("choice ~ ", paste(PredictorVariables, collapse=" + ")))

    multi_fit.1 <- mlogit::mlogit(Formula, data=DM_2, choice="choice", alt.var="altern", chid.var="clust",
                          weights=wgts, shape="long", constPar="offst")

    omegaFullMat <- multi_fit.1$probabilities
    omegaFullMat <- omegaFullMat[,order(as.numeric(colnames(omegaFullMat)))]
    omegaFull <- data.frame(omegaFullMat, subj=ml_dat_trim$subj)
    wtMat <- data.frame(wgts=ml_dat_trim$wgts,subj=ml_dat_trim$subj)
    reWeight <- function(indx, Mat=omegaFull, wt=wtMat){
      subjDat <- Mat[Mat$subj==indx,][,-ncol(Mat)]
      wts <- wt[wt$sub==indx,"wgts"]

      y <- colSums(subjDat*wts)
    }
    omegaFull <- do.call(rbind,lapply(1:nClust, reWeight, Mat=omegaFull, wt=wtMat))
    colnames(omegaFull) <- cellIDs

    uMat <- matrix(0, nrow = nClust, ncol = clustSize)
    for (t in 1:clustSize) {
      uMat[,t] <- rowSums(omegaFull[,substr(cellIDs,t,t)=="1"])
    }
    u <- matrixcalc::vec(uMat)



    # M-step for gamma
    wts <- as.vector(omegaFull)
    omega <- as.vector(t(omegaFull))
    nu2 <- kronecker(diag(R2),rep(1,nClust))
    dfM = cbind.data.frame(z = apply(nu2, 1, function(x) which(x==1)), wts)
    DM <- mlogit::mlogit.data(dfM, shape="wide",sep="",choice="z",alt.levels=1:R2)
    #MstepDat <- data.frame(subj=subj, altern=altern, wts=omega)
    multi_fit.2 <- mlogit::mlogit(z ~ -1|1, weights = wts, data = DM)
    #sv <- round(coef(fit3),2)
    #multi_fit.2 <- myMlogit(parm=sv, multy=MstepDat$wts, xMat=NULL, alt=MstepDat$altern, clust=MstepDat$subj)
    gamma <-coef(multi_fit.2)
    #gamma <- multi_fit.2$newparm
    glpVec <- GMat%*%gamma
    glpMat <- matrix(glpVec, nrow = nClust, ncol = R2-1, byrow = TRUE)

    #Compute the mean for mixing mechanism
    piMat <- getMultip(glpMat)
    piMatFull <- cbind(piMat,1-rowSums(piMat))
    colnames(piMatFull) <- cellIDs # 1 denotes an excess zero
    nu <- piMatFull
    #Setting up the offset
    offst <- log(nu[,1:(ncol(nu))]/nu[,ncol(nu)])[1:nClust,]

    # M step for beta, rho, phi
    if(corstruct == 'indep') {
      glmwts = 1 - u
      if(errDist=='poisson'){
        famfun <- if(nophi) poisson(link=mulink) else quasipoisson(link=mulink)
      } else if(errDist=='binomial'){
        famfun <- if(nophi)  binomial(link=mulink) else quasibinomial(link=mulink)
      }
      glmfit <- suppressWarnings(glm.fit(bmat, y , family = famfun, weights = glmwts))
      beta <- glmfit$coef
      if(!nophi) phi = summary(glmfit)$dispersion
    } else{
      geeits <- 0
      geeConv <- FALSE

      while (geeits < maxits) {
        geeits <- geeits + 1
        blp <- bmat %*% beta
        mu <- getmu(blp, bin = bin, link = mulink)

        dta2 <- data.frame(id = id, time = time, y = y, EY = mu,
                           X = bmat, lp = blp, u = u)
        dtaList2 <- split(dta2, dta2$id)

        # set up the list of data frames and arguments to getScoreHessi
        iList.2 <- lapply(dtaList2, getScoreHessi, corParm = rho, dispParm = phi,
                          nophi = TRUE, errDist = errDist, link = mulink, corstruct = corstruct,
                          weighted = TRUE)

        i2List.2 <- unlist(iList.2, recursive = FALSE)
        udotInd <- grep(".efi", names(i2List.2))
        udotList <- i2List.2[udotInd]
        udot <- Reduce("+", udotList)
        idotInd <- grep(".ii", names(i2List.2))
        idotList <- i2List.2[idotInd]
        idot <- Reduce("+", idotList)

        idotinv <- solve(idot)
        updatevec <- idotinv %*% udot


        maxscore <- max(abs(udot))
        maxchange <- max(abs(updatevec))

        beta <- beta - updatevec[1:dimbeta]
        rho <- rho - updatevec[(dimbeta + 1):(dimbeta + dimrho)]
        phi <- phi - updatevec[dimparmB] * as.numeric(!nophi)
        if(verbose) cat("Beta step ",geeits,"beta=",beta,", rho=",rho,"\n")

        maxparm <- max(c(abs(beta), abs(rho)))
        minparm <- min(c(abs(beta), abs(rho)))

        if (maxparm > BIG || minparm < SMALL || rcond(idot) < SMALL) {
          # not converging
          geeits <- maxits + 1
        }

        if (maxscore < toler && maxchange < toler) {
          geeConv <- TRUE
          geeits <- maxits + 1
        }  #if(!geefail)
      } #End while loop for beta,rho,phi

      if(!geeConv) {
        converge <- FALSE
        stop("S step for beta did not converge")
        #      its = maxits+1
      }
    }
    parm = c(gamma,beta)
    if(corstruct != 'indep'){
      parm <- c(parm,rho)
    }

    if(!nophi){
      parm <- c(parm,phi)
    }

    maxchange <- max(abs(parm-parmold))
    relchange <- vecnorm(parm-parmold)/vecnorm(parm)

    if(relchange < toler){
      converge <- TRUE
    }

  } #End global while loop

  return(list("Gamma"=gamma,"Beta"=beta,"rho"=rho,"phi"=phi,
              "iterations"=its,"converge"=converge))
} #End function

