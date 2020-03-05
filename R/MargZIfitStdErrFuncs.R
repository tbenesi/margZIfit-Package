getScorei <- function(params, index, y, dims, mulink, plink, bmat, gmat, id,
                      time, nClust, clustSize, bin, errDist, corstruct,mixcorstr, nophi){

  # Order of parameters: params <- c(beta, gamma, rho/phi, delta, S)

  # Dimensions
  dimbeta <- dims$dimbeta;  dimgamma <- dims$dimgamma
  dimphi <- dimrho <- dimdelta <- 0
  if(mixcorstr!='indep') dimdelta <- dims$dimdelta
  if(corstruct!='indep') dimrho <- dims$dimrho
  if(!nophi) dimphi <- 1
  dimparm <- dimgamma+dimbeta+dimrho+dimdelta+dimphi

  beta <- params[1:dimbeta]
  gamma <- params[(dimbeta+1):(dimbeta+dimgamma)]
  phi <- 1
  rho <- delta <- 0
  if(corstruct!='indep') rho <- params[(dimbeta+dimgamma+1):(dimbeta+dimgamma+dimrho)]
  if(!nophi) phi <- params[(dimbeta+dimgamma+dimrho+1)]
  if(mixcorstr!='indep') delta<-params[(dimbeta+dimgamma+dimrho+dimphi+1):(dimbeta+dimgamma+dimrho+dimphi+dimdelta)]


  # Compute the weights
  bmat_i <- bmat[id==index,]
  blp_i <- as.matrix(bmat_i)%*%beta
  gmat_i <- as.matrix(gmat[id==index,])
  glp_i <- gmat_i%*%gamma
  pi <- getmixp(glp_i, link = plink)
  Bin <- bin[id==index]
  mu2 <- getmu(blp_i, bin = Bin, link = mulink)
  yi <- y[id==index]
  ni <- length(yi)
  ui <- rep(0, ni)
  fofyij <- 1 + ((1-pi)/pi)*getf2(yi, mu = mu2, bin = Bin, errDist = errDist)
  ui[yi==0] <- 1/fofyij[yi==0]

  id_i <- id[id==index]
  time_i <- time[id==index]
  dta <- data.frame(id = id_i, time = time_i, y = ui, EY = pi, X = gmat_i, lp = glp_i)

  # set up the list of data frames and arguments to getScoreHessi
  i2List <- getScoreHessi(dta, corParm = delta, dispParm = 1, nophi = TRUE, errDist = "binomial",
                          link = plink, corstruct = mixcorstr)
  gamma_part <- i2List$gamma_part
  delta_part <- i2List$delta_part

  dta2 <- data.frame(id = id_i, time = time_i, y = yi, EY = mu2, X = bmat_i, lp = blp_i, u = ui)
  iList.2 <- getScoreHessi(dta2, corParm = rho, dispParm = phi, nophi = TRUE,
                           errDist = errDist, link = mulink, corstruct = corstruct, weighted = TRUE)
  beta_part <- iList.2$beta_part
  rho_part <- iList.2$rho_part

  Ui <- c(beta_part, gamma_part)
  if(corstruct!='indep') Ui <- c(Ui, rho_part)
  if(mixcorstr!='indep') Ui <- c(Ui, delta_part)
  return(Ui)
}



getMultiScorei <- function(params, Index, y, dims, mulink, plink, bmat, gmat, id,
                           time, nClust, clustSize, bin, errDist, corstruct,nophi){

  # Order of parameters: params <- c(beta, gamma, rho/phi, delta, /S)

  # Dimensions
  dimbeta <- dims$dimbeta;  dimgamma <- dims$dimgamma
  dimphi <- dimrho <- 0
  if(corstruct!='indep') dimrho <- dims$dimrho
  if(!nophi) dimphi <- 1
  dimparm <- dimgamma+dimbeta+dimrho+dimphi

  # Paramerters
  beta <- params[1:dimbeta];  gamma <- params[(dimbeta+1):(dimbeta+dimgamma)]
  phi <- 1;  rho <- 0
  if(corstruct!='indep') rho <- params[(dimbeta+dimgamma+1):(dimbeta+dimgamma+dimrho)]
  if(!nophi) phi <- params[(dimbeta+dimgamma+dimrho+1)]

  N <- nClust*clustSize;  R <- 2^clustSize
  #Compute the conditional mixing probs
  y0 <- as.numeric(y==0)
  y0Mat <- matrix(y0, nrow = nClust, ncol = clustSize, byrow = FALSE)
  y0Tab <- ftable(as.data.frame(y0Mat), col.vars = 1:clustSize)
  Noclass <- which(y0Tab==0)
  len_Noclass <- length(Noclass); R2 <- R-len_Noclass
  if(len_Noclass==0) Noclass <- 1000*R
  y0Tab <- y0Tab[,-Noclass]


  cellids <- R.utils::intToBin(0:(R-1))[-Noclass]

  ##Compute the E-step for multinomial mixing
  #Compute the linear predictor for the non-degenerate distribution
  blp = bmat%*%beta
  bin = bin
  #Compute the mean for the non-degenerate component
  mu <- getmu(blp, bin = bin, link = mulink)

  #initialize the estimates of the conditional mixing prob
  GMatList <- lapply(1:nClust, function(x) kronecker(diag(R2-1),gmat[x,]))
  GMat <- Reduce(rbind,GMatList, accumulate = FALSE)


  #GMat <- getMultiGMAT(gmat,nClust,Index)
  idvec <- kronecker(unique(id),rep(1,R2-1))
  glpVec <- GMat%*%gamma
  glpMat <- matrix(glpVec, nrow = nClust, ncol = R2-1, byrow = TRUE)

  #Compute the mean for mixing mechanism
  piMat <- getMultip(glpMat)
  piMatFull <- cbind(piMat,1-rowSums(piMat))
  colnames(piMatFull) <- cellids # 1 denotes an excess zero
  nu <- piMatFull

  f2ofy <- getf2(y, mu, bin = bin, errDist = errDist)
  f2ofyMat <- matrix(f2ofy, nrow = nClust, ncol = clustSize, byrow = FALSE)

  EstepTerms <- matrix(0, nrow = nClust, ncol = R2)
  cellids <- R.utils::intToBin(0:(R-1))[-Noclass]
  for(j in 1:R2){
    indices <- unlist(gregexpr("0",cellids[j]))
    if(all(indices==-1)) indices <- -(1:clustSize)
    EstepTerms[,j] <- nu[,j,drop=FALSE]*rowProds(f2ofyMat[,indices,drop=FALSE])*
      rowProds(y0Mat[,-indices,drop=FALSE])
  }
  omegaFull <- EstepTerms/rowSums(EstepTerms)

  uMat <- matrix(0, nrow = nClust, ncol = clustSize)
  for (t in 1:clustSize) {
    uMat[,t] <- rowSums(omegaFull[,substr(cellids,t,t)=="1"])
  }
  u <- matrixcalc::vec(uMat)

  yi <- y[id==Index]; id_i <- id[id==Index]; ui <- u[id==Index]
  mui <- mu[id==Index]; time_i <- time[id==Index]
  bmat_i <- bmat[id==Index,]; blp_i <- blp[id==Index]

  dta2 <- data.frame(id = id_i, time = time_i, y = yi, EY = mui,
                     X = bmat_i, lp = blp_i, u = ui)
  iList.2 <- getScoreHessi(dta2, corParm = rho, dispParm = phi,
                           nophi = TRUE, errDist = errDist, link = mulink, corstruct = corstruct,
                           weighted = TRUE)

  beta_part <- iList.2$beta_part
  rho_part <- iList.2$rho_part

  glpVec <- GMat%*%gamma
  YVec = matrixcalc::vec(t(omegaFull[,-R2]))
  glpMat <- matrix(glpVec, nrow = nClust, ncol = R2-1, byrow = TRUE)
  piMat <- getMultip(glpMat)
  piVec <- matrixcalc::vec(t(piMat))
  idVec <- kronecker(unique(id),rep(1,R2-1))

  idVec_i <- idVec[idVec==Index]; glpVec_i = glpVec[idVec==Index]
  YVec_i <- YVec[idVec==Index]; piVec_i <- piVec[idVec==Index] ;
  GMat_i <-as.matrix(GMat[idVec==Index, ])
  dta <- data.frame(id=idVec_i, glpVec=glpVec_i, y=YVec_i, EYvec=piVec_i, X = GMat_i)

  # set up the list of data frames and arguments to getScoreHessi
  iList <- MultinomFuncs(dta)

  gamma_part <- iList$efdot

  Ui <- c(beta_part, gamma_part)
  if(corstruct!='indep') Ui <- c(Ui, rho_part)
  return(Ui)
} #End function


getMultiMBScorei <- function(params, Index, y, omega=omega, u=u, dims, mulink, plink, bmat,
                             gmat, id, time, nClust, clustSize, bin, errDist, corstruct, nophi){

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
  #dimgamma <- dimgamma - len_Noclass
  #gamma <- gamma[-Noclass]
  cellIDs <- R.utils::intToBin(0:(R-1))[-Noclass] # 1 denotes an excess zero
  R2 <- R-len_Noclass

  GMatList <- lapply(1:nClust, function(x) kronecker(diag(R2-1),gmat[x,]))
  GMat <- Reduce(rbind,GMatList, accumulate = FALSE)


  #   # Compute the linear predictor for the non-degenerate distribution
  blp <- bmat%*%beta
  # Compute the mean for the non-degenerate component
  mu <- getmu(blp, bin = bin, link = mulink)
  omegaFull <- omega; u <- u

  yi <- y[id==Index]; id_i <- id[id==Index]; ui <- u[id==Index]
  mui <- mu[id==Index]; time_i <- time[id==Index]
  bmat_i <- bmat[id==Index,]; blp_i <- blp[id==Index]

  dta2 <- data.frame(id = id_i, time = time_i, y = yi, EY = mui,
                     X = bmat_i, lp = blp_i, u = ui)
  # set up the list of data frames and arguments to getScoreHessi
  iList.2 <- getScoreHessi(dta2, corParm = rho, dispParm = phi,
                           nophi = nophi, errDist = errDist, link = mulink, corstruct = corstruct,
                           weighted = TRUE)

  beta_part <- iList.2$beta_part
  rho_part <- iList.2$rho_part

  glpVec <- GMat%*%gamma
  YVec = matrixcalc::vec(t(omegaFull[,-R2]))
  glpMat <- matrix(glpVec, nrow = nClust, ncol = R2-1, byrow = TRUE)
  piMat <- getMultip(glpMat)
  piVec <- matrixcalc::vec(t(piMat))
  idVec <- kronecker(unique(id),rep(1,R2-1))

  idVec_i <- idVec[idVec==Index]; glpVec_i = glpVec[idVec==Index]
  YVec_i <- YVec[idVec==Index]; piVec_i <- piVec[idVec==Index] ;
  GMat_i <-as.matrix(GMat[idVec==Index, ])
  dta <- data.frame(id=idVec_i, glpVec=glpVec_i, y=YVec_i, EYvec=piVec_i, X = GMat_i)

  # set up the list of data frames and arguments to getScoreHessi
  iList <- MultinomFuncs(dta)

  gamma_part <- iList$efdot
  Ui <- c(beta_part, gamma_part)
  if(corstruct!='indep') Ui <- c(Ui, rho_part)

  return(Ui)
} #End function


getEstep <- function(params, y, dims, mulink, plink, bmat, gmat, id, time, nClust, clustSize,
                     bin, errDist, corstruct, nophi){
  dimbeta <- dims$dimbeta; dimgamma <- dims$dimgamma
  N = nClust*clustSize
  R <- 2^clustSize

  dimrho <- dimphi <- 0
  if(corstruct!='indep') dimrho <- dims$dimrho
  if(!nophi) dimphi <- 1

  # The response matrix
  Ymat <- matrix(y, ncol = clustSize , byrow = FALSE)
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
  cellids <- R.utils::intToBin(0:(R-1))[-Noclass] # 1 denotes an excess zero
  R2 <- R-len_Noclass
  y0Tab <- y0Tab[-Noclass]
  nuFullInit <- matrix(c(y0Tab/sum(y0Tab)),nrow = nClust, ncol = R2, byrow = TRUE)
  colnames(nuFullInit) <- cellids
  nuInit <- nuFullInit[,1:(R2-1)]
  gamma <- gamma[-Noclass]

  # set up initial values for the omega_ir's. These values are indicators for whether the
  # rth pattern of zero valued responses has occurred in the ith cluster.
  # so omega_ir is not constant over i
  yi0Tab <- ftable(data.frame(y0Mat, sub = 1:nClust),col.vars = 1:clustSize, row.vars = "sub")[,-Noclass]
  omegaFullInit <- as.matrix(yi0Tab)
  omegaInit <- omegaFullInit[,1:(R2-1)]

  parm <- c(gamma,beta)
  if(corstruct!='indep') parm <- c(parm,rho)
  if(!nophi) parm <- c(parm,phi)


  # Compute the linear predictor for the non-degenerate distribution
  blp <- bmat%*%beta
  # Compute the mean for the non-degenerate component
  mu <- getmu(blp, bin = bin, link = mulink)

  # Estimates of the conditional mixing probs
  nu <- nuFullInit #set nu to initial estimates
  offst <- log(nuFullInit[,1:(ncol(nuFullInit))]/nuFullInit[,ncol(nuFullInit)])
  omega = omegaFullInit


  # Setup the means to compute the covariates
  muMat <- matrix(mu, nrow = nClust , ncol = clustSize)
  Mulist <- split(muMat, 1:nrow(muMat))
  all_pos <- lapply(strsplit(cellids,''),function(a) 1-as.numeric(a))
  all_posList <- rep(all_pos, times=nClust)
  Vechall_pos <- lapply(all_posList, function(x) vech(x%*%t(x)))
  Ylist <- split(Ymat, 1:nrow(Ymat))
  riList <- purrr::map2(Ylist, Mulist, function(x,y) x-y)
  #cross_riList <- lapply(riList, prod)
  riList_rep <- rep(riList,each=R2)
  #cross_riList_rep <- rep(cross_riList,each=R2)
  #riList2 <- map2(riList_rep,cross_riList_rep,c)
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


  #Compute the conditional mixing probs (E-Step:)
  f2ofy <- getf2(y, mu, bin = bin, errDist = errDist)
  f2ofyMat <- matrix(f2ofy, nrow = nClust, ncol = clustSize, byrow = FALSE)
  EstepTerms <- matrix(0, nrow = nClust, ncol = R2)
  cellids <- R.utils::intToBin(0:(R-1))[-Noclass]
  for(j in 1:R2){
    indices <- unlist(gregexpr("0",cellids[j]))
    if(all(indices==-1)) indices <- -(1:clustSize)
    EstepTerms[,j] <- nu[,j,drop=FALSE]*rowProds(f2ofyMat[,indices,drop=FALSE])*
      rowProds(y0Mat[,-indices,drop=FALSE])
  }
  omega <- EstepTerms/rowSums(EstepTerms)
  wgts <- matrixcalc::vec(t(omega))


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

  tryCatch({multi_fit.1 <- mlogit::mlogit(Formula, data=DM_2, choice="choice", alt.var="altern", chid.var="clust",
                                  weights=wgts, shape="long", constPar="offst")},
           error=function(e){multi_fit.1 <<- NULL})

  if(is.null(multi_fit.1)){
    xMat <- as.matrix(cbind(ml_dat_trim[, colnames(ml_dat_trim)[grepl('X',colnames(ml_dat_trim))]],
                            ml_dat_trim$offst))
    sv = rep(1,nrow(xMat))
    multi_fit.1 <- myMlogit(parm=sv, multy=ml_dat_trim$wgts, alt=ml_dat_trim$altern, xMat=xMat,
                            clust=ml_dat_trim$subj, yMat=Ymat, muMat=muMat)
    omegaFull <- multi_fit.1$piMat
  }else{
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
  }
  colnames(omegaFull) <- cellids
  uMat <- matrix(0, nrow = nClust, ncol = clustSize)
  for (t in 1:clustSize) {
    uMat[,t] <- rowSums(omegaFull[,substr(cellids,t,t)=="1"])
  }
  u <- matrixcalc::vec(uMat)

  return(list(omega=omegaFull,u=u))
}


