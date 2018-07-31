#' r_mediation : function to simulate DNA methylation data for mediation analyze
#'
#' @param n	: number of individuals
#' @param p	: number of cpg variables
#' @param K	: number of latent factors
#' @param freq : (vector) mean methylation values (if NULL, set randomly)
#' @param prop.causal : proportion of causal variables (probes/loci)
#' @param prop.causal.x : proportion of causal cpg M -> x
#' @param prop.causal.y : proportion of causal cpg M -> y
#' @param prop.causal.ylx : proportion of causal y in causal x
#' @param prop.variance.y : proportion of phenotypic variance explained by latent structure (intensity of confounding)
#' @param prop.variance.x : proportion of exposure variance explained by latent structure (intensity of confounding)
#' @param rho : correlation outcome/exposure (direct effect)
#' @param sigma :	standard deviation of residual errors
#' @param sd.B : standard deviation for effect sizes (B: M->Y)
#' @param mean.B :	(vector) mean of effect sizes
#' @param sd.A :	standard deviation for effect sizes (A: M->X)
#' @param mean.A :	(vector) mean of effect sizes
#' @param sd.U : (vector) standard deviations for factors
#' @param sd.V : standard deviations for loadings
#'
#' @return M : matrix of methylation beta values
#' @return Y : phenotype/health outcome
#' @return B : effect sizes phenotype/health outcome
#' @return X : exposure
#' @return A : effect sizes exposure
#' @return mediators : set of true mediators
#' @return causal.x : set of CpGs associated with the exposure
#' @return causal.y : set of CpGs associated with the outcome
#' @return U : simulated confounders
#' @return V : loadings of coufounders
#' @return freq : mean methylation values
#' @return controls : true control gene (NOT USE for simulation study)
#'
#' @details
#'
#' This function is used to simulate datasets for analysis of mediations.
#' The simulation model is based on linear relationships.
#' First, it construct a covariance matrix for X, Y and U using the parameter rho
#' (direct effect or correlation between X and Y) and propvar
#' (intensity of the confounders or correlation between Y and U).
#' Then this matrix is used to simulate via normal laws X, Y and U.
#' Thereafter, the effect sizes of X (A), Y (B) and U (V) are calculated
#' using mean parameters of effect sizes (meanA and meanB) and standard deviations (sdA, sdB and sdV).
#' Note that the effect sizes of X and Y are calculated only for causal mediators with X and/or Y.
#' For non-causal mediators, the effect sizes is 0.
#' On the other hand, a residual error matrix is calculated via the sigma (Z) parameter.
#' To finish the methylation matrix is calculated thanks to the formula : M = V*U + A*X + B*Y + Z
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' @export
r_mediation <- function(n,
                         p,
                         K,
                         freq = NULL,
                         prop.causal.x = 0.010,
                         prop.causal.y = 0.010,
                         prop.causal.ylx = 0.5,
                         prop.variance.y = 0.6,
                         prop.variance.x = 0.2,
                         rho = 0.2,
                         sigma = 0.2,
                         sd.A = 1.0,
                         mean.A = 3.0,
                         sd.B = 1.0,
                         mean.B = 5.0,
                         sd.U = 1.0,
                         sd.V = 1.0)
{
  causal.x <- sample.int(p, prop.causal.x * p)
  causal.ylx <- sample(causal.x , prop.causal.ylx*length(causal.x))
  if (prop.causal.y * p < prop.causal.ylx*length(causal.x)) {
    stop("# causal y < # mediators")
  }
  else {
    causal.y <- c(causal.ylx, sample.int(p, prop.causal.y * p - prop.causal.ylx*length(causal.x)) )
  }
  x.nb = length(causal.x)
  y.nb = length(causal.y)

  if (is.null(freq)) freq <- runif(n = p,min =  0.2,max =  0.8) # mean of methylation for each site

  if (prop.variance.y + rho^2 > 1) stop("prop.variance.y + rho^2 > 1")
  if (prop.variance.x + rho^2 > 1) stop("prop.variance.x + rho^2 > 1")

  #cs <- runif(K, min = -1, max = 1)
  #theta.y <- sqrt( prop.variance.y /sum((cs/sd.U)^2) )
  #theta.x <- sqrt( prop.variance.x /sum((cs/sd.U)^2) )

  cs.y <- runif(K, min = -1, max = 1)
  cs.x <- runif(K, min = -1, max = 1)
  theta.y <- sqrt( prop.variance.y /sum((cs.y/sd.U)^2) )
  theta.x <- sqrt( prop.variance.x /sum((cs.x/sd.U)^2) )

  # constructing the covariance matrix
  Sigma <- diag(x = sd.U^2, nrow = K, ncol = K)

  Sigma <- rbind(Sigma, matrix(cs.y*theta.y, nrow = 1))
  Sigma <- rbind(Sigma, matrix(cs.x*theta.x, nrow = 1))

  Sigma <- cbind(Sigma, matrix(c(cs.y*theta.y, 1, rho), ncol = 1))
  Sigma <- cbind(Sigma, matrix(c(cs.x*theta.x, rho, 1), ncol = 1))

  UYX <- MASS::mvrnorm(n, mu = rep(0, K + 2), Sigma = Sigma)
  U <- UYX[, 1:K, drop = FALSE]   # confounders
  Y <- UYX[, K + 1, drop = FALSE] # outcome
  X <- UYX[, K + 2, drop = FALSE] # exposure

  V <- MASS::mvrnorm(p, mu = rep(0, K), Sigma = sd.V^2 * diag(K))

  A <- matrix(0, p, 1)
  A[causal.x, 1] <- rnorm(x.nb, mean.A, sd.A)

  B <- matrix(0, p, 1)
  B[causal.y, 1] <- rnorm(y.nb, mean.B, sd.B)

  Epsilon <- apply(matrix(rep(0,p),nrow = 1), 2, function(x) rnorm(n,x,sigma))

  Z = U %*% t(V) + X %*% t(A) + Y %*% t(B) + Epsilon

  M = matrix(rep(qnorm(freq),n), nrow = n, byrow = T) + Z

  M = pnorm(M)

  return(list(M = M,
              Y = Y,
              B = B,
              X = X,
              A = A,
              mediators = sort(causal.ylx),
              causal.x = sort(causal.x),
              causal.y = sort(causal.y),
              U = U,
              V = V,
              freq = freq,
              Sigma = Sigma,
              controls = !(1:p %in% unique(sort(c(sort(causal.x),sort(causal.y)))))))
}



#' wrapper of mediation package (for one mediator)
#'
#' @param X exposure (continuous)
#' @param Y outcome (continuous or factor)
#' @param mediator one CpG site
#' @param CONF estimated and known confounders
#' @param sims number of simulations
#' @param ... other options for the mediation model
#'
#' @return a pValue for ACME (ACME.p.med)
#' @return a effect for ACME (ACME.eff)
#'
#' @details
#'
#' This function makes it possible to carry out an independent mediation analysis between an exposure (X),
#'  a CpG (M) and an outcome (Y). For continuous X and Y it carries out a conbinaison of 2 linear models.
#'   For a continuous X and a factorial Y it is a conbinaison of a linear model and a GLM.
#'   With these conbinations, it calculates an ACME (Average Causal Mediation Effects) and an associated value.
#'
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # test if the first CpG mediate a effect
#' mod <- wrapMed(simu$X, simu$Y, simu$M[,1])
#' @export
wrapMed <- function(X,Y,mediator,CONF = NULL,sims=100) {

  if (is.null(CONF) == T) { # if no confouders

    da.med <- data.frame(m = mediator, x = X, y = Y)

    med1 <- lm(m ~ x,data = da.med)

    if (is.factor(Y)) {
      med2 <- glm(y ~ m + x, data = da.med,family = binomial("probit"))
    }
    else {
      med2 <- lm(y ~ m + x, data = da.med)
    }

    out1 <- mediation::mediate(med1,med2,treat = "x",mediator = "m",sims = sims)

  }
  else {
    da.med <- data.frame(m = mediator, x = X, y = Y, conf = CONF)

    med1 <- lm(m ~ .,data = da.med[, -3])

    if (is.factor(Y) == TRUE) {
      med2 <- glm(y ~ ., data = da.med,family = binomial("probit"))
    }
    else {
      med2 <- lm(y ~ ., data = da.med)
    }

    out1 <- mediation::mediate(med1,med2,treat = "x",mediator = "m",sims = sims)

  }

  ACME.e <- out1$d.avg
  ACME.p <- out1$d.avg.p

  # other test

  sim <- c(as.vector(out1$d0.sims),as.vector(out1$d1.sims))

  return(list(ACME.eff = ACME.e, ACME.p.med = ACME.p, sim = sim))
}


#' wrapGLM : naive glm function for mediation EWAS and independent mediation analysis (no estimated confondeurs)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param family "binomial" or "gaussian" : family function for GLM
#' @param sims Number of MC simulations in the mediation model
#'
#' @return pval, qval, score, standard errors and effect sizes for each CpG
#' 
#' @details
#' 
#' This function uses a generalized linear model "glm" to perform association tests of the mediation EWAS. 
#' Note that no confusers are estimated in this function.
#' 
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis
#' mod <- wrapGLM(simu$X, simu$Y, mediation = T)
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#' @export
wrapGLM <- function(X,Y,M,conf.known=NULL,mediation=F,top=20,family="binomial",sims = 100) {

  if (is.null(conf.known) == T) {
    obj <- data.frame(X, Y)
  }
  else {
    obj <- data.frame(X, Y,conf.known)
  }

  if (family == "binomial") {
    coeff <- t(apply(M,
                     2,
                     function(x) summary(glm(x ~ .,
                                             data = obj,
                                             family = binomial(link = "probit")))$coeff[2:(1+ncol(obj)), 1:4]))
  }
  else {
    coeff <- t(apply(M,
                     2,
                     function(x) summary(glm(x ~ .,
                                             data = obj))$coeff[2:(1+ncol(obj)), 1:4]))

  }

  es <- coeff[, 1 : ncol(obj)]
  se <- coeff[, (1+ncol(obj)) : (2*ncol(obj))]
  sc <- coeff[, (1+2*ncol(obj)) : (3*ncol(obj))]
  pv <- coeff[, (1+3*ncol(obj)) : (4*ncol(obj))]


  # pValue unique
  pval <- apply(pv[, 1:2], 1, max)^2
  # qValue by FDR
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F)$qval

  # mediation
  if (mediation == T) {
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = conf.known,sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)

  }
  else {
    media <- NULL
  }

  return(list(pval   = pval,
              qval   = qval,
              score  = sc[, 1:2],
              es.bac = es,
              se.bac = se,
              media  = media))
}

#' wrapCATE : Mediation EWAS and independent mediation analysis (coufounders estimated via CATE method)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param K number of coufounders to estimate
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param sims Number of MC simulations in the mediation model
#' @param ... other options of the CATE model
#'
#'
#' @return pval, qval, score for each CpG (note: the pValues are calibrated)
#' @return conf : coufounders estimated via CATE method
#' @return media : if mediation = T, result of wrapMed() function
#' @return mod : whole CATE model
#'
#' @details 
#' 
#' This function estimates the confusers via the CATE model. 
#' In addition, the CATE model uses robust linear regressions "rlm" to perform EWAS associations testing of mediation.
#' 
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis
#' mod <- wrapCATE(simu$X, simu$Y, 5, mediation = T)
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#' @export
wrapCATE <- function(X,Y,M,K,conf.known = NULL,mediation = T,top = 50,sims = 100, ...) {

  if (is.null(conf.known) == T) {
    mod <- cate::cate.fit(X.primary = cbind(X, Y), X.nuis = rep(1,nrow(cbind(X, Y))), Y = M, r = K)
  }
  else {
    mod <- cate::cate.fit(X.primary = cbind(X, Y), X.nuis = cbind(1,conf.known), Y = M, r = K)
  }

  # pValue unique
  pval <- apply(mod$beta.p.value[, 1:2], 1, max)^2
  # qValue by FDR
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F,...)$qval
  score.x <- mod$beta.t[,1]
  score.y <- mod$beta.t[,2]

  # confounders
  conf <- mod$Z

  # mediation
  if (mediation == T) {
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = cbind(conf,conf.known),sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)
  }
  else {
    media <- NULL
  }

  return(list(pval  = pval,
              qval  = qval,
              conf  = conf,
              score = cbind(score.x,score.y),
              mod   = mod,
              media = media))
}

#' wrapLFMM : Mediation EWAS and independent mediation analysis (coufounders estimated via LFMM methods)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param K number of coufounders to estimate
#' @param meth LFMM combines two methods to estimate coufondeurs: "ridge" and "lasso"
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param sims Number of MC simulations in the mediation model
#' @param ... other options of the LFMMs models
#'
#' @return pval, qval, score for each CpG (note: the pValues are calibrated)
#' @return conf : coufounders estimated via LFMM methods
#' @return media : if mediation = T, result of wrapMed() function
#' @return mod : whole LFMM model
#' 
#' @details 
#' 
#' This function estimates the confusers via the LFMM model. There are two sub-models of LFMM to estimate the confondeurs, 
#' "lasso" and "ridge".  
#' In addition, the LFMM model uses linear regression "lm" to perform associations tests of mediation EWAS.
#'
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis (LFMM ridge)
#' mod <- wrapLFMM(simu$X, simu$Y, 5, mediation = T,meth = "ridge")
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#' @export
wrapLFMM <- function(X,Y,M,K,meth=c("ridge","lasso"),conf.known = NULL,mediation = T,top = 50,sims = 100,...) {

  if (is.null(conf.known) == T) {
    Xnew <- cbind(X,Y)
  }

  else {
    Xnew <- cbind(X,Y,conf.known)
  }

  if (meth == "ridge") {
    # LFMM ridge
    mod <- lfmm::lfmm_ridge(Y = M, X = Xnew, K = K,...)
  }
  else {
    # LFMM lasso
    mod <- lfmm::lfmm_lasso(Y = M, X = Xnew, K = K,...)
  }

  mod1 <- lfmm::lfmm_test(Y = M, X = Xnew, lfmm = mod)
  # pValue unique
  pval <- apply(mod1$calibrated.pvalue[, 1:2], 1, max)^2 # pval unique
  # qValue by FDR
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F)$qval
  score.x <- mod1$score[,1]
  score.y <- mod1$score[,2]
  # confounders
  conf <- mod$U

  # mediation
  if (mediation == T) {
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = cbind(conf,conf.known),sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)
  }
  else {
    media <- NULL
  }

  return(list(pval  = pval,
              qval  = qval,
              conf  = conf,
              score = cbind(score.x,score.y),
              mod   = list(mod = mod,test = mod1),
              media = media))
}

#' wrapDSVA : Mediation EWAS and independent mediation analysis (coufounders estimated via dSVA method)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param K number of coufounders to estimate
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param sims Number of MC simulations in the mediation model
#'
#' @return pval, qval, score for each CpG
#' @return conf : coufounders estimated via dSVA method
#' @return media : if mediation = T, result of wrapMed() function
#' @return mod : whole dSVA model
#' 
#' @details 
#' 
#' This function estimates the confusers via the dSVA model. 
#' In addition, the dSVA model use linear regression "lm" to perform associations tests of mediation EWAS.
#'
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis
#' mod <- wrapDSVA(simu$X, simu$Y, 5, mediation = T)
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#' @export
wrapDSVA <- function(X,Y,M,K,conf.known = NULL,mediation = T,top = 50, sims =100) {
  if (is.null(conf.known) == T) {
    Xnew <- cbind(X,Y)
  }

  else {
    Xnew <- cbind(X,Y,conf.known)
  }
  # dSVA model
  mod <- dSVA::dSVA(Y = M,X = Xnew,ncomp = K)
  # pValue unique
  pval <- apply(t(mod$Pvalue[1:2,]), 1, max)^2
  # qValue by FDR
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F)$qval
  score.x <- mod$Bhat[1,]/mod$BhatSE[1,]
  score.y <- mod$Bhat[2,]/mod$BhatSE[2,]
  # confounders
  conf <- mod$Z

  # mediation
  if (mediation == T) {
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = cbind(conf,conf.known), sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)
  }
  else {
    media <- NULL
  }

  return(list(pval  = pval,
              qval  = qval,
              conf  = conf,
              score = cbind(score.x,score.y),
              mod   = mod,
              media = media))
}


#' wrapSVA : Mediation EWAS and independent mediation analysis (coufounders estimated via SVA method)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param K number of coufounders to estimate
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param controls NOT USE
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param sims Number of MC simulations in the mediation model
#' @param ... other options of the SVA model
#'
#' @return pval, qval, score for each CpG
#' @return conf : coufounders estimated via SVA method
#' @return media : if mediation = T, result of wrapMed() function
#' @return mod : whole SVA model
#' 
#' #' @details 
#' 
#' This function estimates the confusers via the SVA model. 
#' In addition, the SVA model use linear regression "lm" to perform associations tests of mediation EWAS.
#'
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis
#' mod <- wrapSVA(simu$X, simu$Y, 5, mediation = T)
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#' @export
wrapSVA <- function(X,Y,M,K,controls = NULL,conf.known = NULL,mediation = T,top = 50 ,sims = 100,...) {
  # SVA model
  if (is.null(conf.known) == T) {
    simuXY <- data.frame(simuX = X, simuY = Y)
    mod1   <- model.matrix(object = ~ simuX + simuY, data = simuXY)
    mod0   <- model.matrix(object = ~ 1, data = simuXY)
  }

  else {
    simuXY <- data.frame(simuX = X, simuY = Y, conf.known)
    conf.known <- data.frame(conf.known)
    mod1   <- model.matrix(object = ~ ., data = simuXY)
    mod0   <- model.matrix(object = ~ ., data = conf.known)
  }

  mod    <- sva::sva(t(M), mod = mod1, mod0 = mod0, n.sv = K,...)
  modSv  <- cbind(mod1,mod$sv)
  mod0Sv <- cbind(mod0,mod$sv)
  # pValue

  x.mod0Sv <- cbind(mod0Sv, modSv[,2])
  y.mod0Sv <- cbind(mod0Sv, modSv[,3])

  px <- sva::f.pvalue(t(M),modSv, y.mod0Sv)
  py <- sva::f.pvalue(t(M),modSv, x.mod0Sv)

  pval <- apply(cbind(px, py), 1, max)^2

  # qValue by FDR
  qval  <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F)$qval

  # score
  sx <- sva::fstats(t(M),modSv, y.mod0Sv)
  sy <- sva::fstats(t(M),modSv, x.mod0Sv)

  sx <- as.vector(sx)
  sy <- as.vector(sy)

  # confounders
  conf <- mod$sv

  # mediation
  if (mediation == T) {
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = cbind(conf,conf.known), sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)
  }
  else {
    media <- NULL
  }

  return(list(pval  = pval,
              qval  = qval,
              conf  = conf,
              score = cbind(sx, sy),
              mod   = mod,
              media = media))
}


#' wrapRUV : Mediation EWAS and independent mediation analysis (coufounders estimated via RUV methods)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param K number of coufounders to estimate
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param ctl RUV methods need a control gene vector. if it is not known, it will be estimated by naive glms
#' @param r.ctl alpha risk for the estimation of the control gene vector
#' @param meth RUV combines two methods to estimate coufondeurs: "RUV2" and "RUV4"
#' @param GLMwrap result of wrapGLM to estimate the control gene
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param sims Number of MC simulations in the mediation model
#' @param ... other options of the RUVs models
#'
#' @return pval, qval, score for each CpG
#' @return conf : coufounders estimated via RUV methods
#' @return media : if mediation = T, result of wrapMed() function
#' @return mod : whole RUV model
#' 
#' #' @details 
#' 
#' This function estimates the confusers via the RUV model. There are two sub-models of RUV to estimate the confondeurs, 
#' "RUV2" and "RUV4". Note that the RUV model uses control genes to estimate confounders. 
#' If you do not know it, it can be estimated via the wrapGLM function.
#' In addition, the RUV model use linear regression "lm" to perform associations tests of mediation EWAS.
#'
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis
#' mod <- wrapRUV(simu$X, simu$Y, 5, mediation = T, meth ="RUV2)
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#'
#' # or with a previous glm
#' mod1 <- wrapGLM(simu$X, simu$Y, mediation = F)
#' mod2 <- wrapRUV(simu$X, simu$Y, 5, mediation = T, meth ="RUV2, GLMwrap = mod1)
#' @export
wrapRUV <- function(X,Y,M,K,ctl=NULL,r.ctl=0.2,meth=c("RUV2","RUV4"),conf.known=NULL,mediation=T,top=50,
                    GLMwrap = NULL, sims = 100,...) {
  # search gene control
  if (is.null(ctl) == T & is.null(GLMwrap) == T) {
    pv.ctl <- wrapGLM(X,Y,M)$pval
    ctl <- ifelse(pv.ctl < r.ctl,F,T)
  }
  else {
    if (is.null(GLMwrap) == F){
      pv.ctl <- GLMwrap$pval
      ctl <- ifelse(pv.ctl < r.ctl,F,T)
    }
  }

  if (meth == "RUV2") {
    # RUV2 model
    mod <- ruv::RUV2(Y = M, X = cbind(X,Y), ctl = ctl, k = K, Z = conf.known,...)
  }
  else {
    # RUV4 model
    mod <- ruv::RUV4(Y = M, X = cbind(X,Y), ctl = ctl, k = K, Z = conf.known,...)
  }

  # pValue unique
  pval <- apply(t(mod$p), 1, max)^2
  # qValue by FDR
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F)$qval
  score.x <- t(mod$t)[,1]
  score.y <- t(mod$t)[,2]
  # Confounders
  conf <- mod$W

  # mediation
  if (mediation == T) {
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = cbind(conf,conf.known), sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)
  }
  else {
    media <- NULL
  }

  return(list(pval  = pval,
              qval  = qval,
              conf  = conf,
              score = cbind(score.x,score.y),
              mod   = mod,
              media = media))
}

#' wrapRFE : Mediation EWAS and independent mediation analysis (coufounders estimated via RefFreeEWAS methods)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param K number of coufounders to estimate
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param nbboot number of bootstrap to calculate pValue
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param sims Number of MC simulations in the mediation model
#' @param ... other options of the RefFreeEWAS model
#'
#' @return pval, qval, score for each CpG
#' @return conf : coufounders estimated via RefFreeEWAS methods
#' @return media : if mediation = T, result of wrapMed() function
#' @return mod : whole RefFreeEWAS model
#' 
#' #' @details 
#' 
#' This function estimates the confusers via the RefFreeEWAS model. 
#' In addition, the RefFreeEWAS model use a bootstrap method to 
#' calculate the pValues of the associations tests of mediation EWAS.
#' Note that it is preferable to use confondeurs estimate by RefFreeEWAS in another test for association (glm for example) 
#' than using pValues calculated by the bootstrap method RefFreeEWAS.
#'
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis
#' mod <- wrapRFE(simu$X, simu$Y, 5, mediation = T)
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#' @export
wrapRFE <- function(X,Y,M,K,nbboot=50,conf.known = NULL,mediation = T,top = 50, sims = 100,...) {
  # RefFreeEwas model
  mod <- RefFreeEWAS::RefFreeEwasModel(t(M),cbind(1, X, Y, conf.known),K = K)
  # bootstrap method
  testBoot <- RefFreeEWAS::BootRefFreeEwasModel(mod,nbboot)
  # score calculation
  score.x <- apply(testBoot, 1:3, mean)[,2,1]/apply(testBoot, 1:3, sd)[,2,1]
  score.y <- apply(testBoot, 1:3, mean)[,3,1]/apply(testBoot, 1:3, sd)[,3,1]
  # pValues calculation
  pvalu.x <- pchisq(score.x^2,df = 1, low = F)
  pvalu.y <- pchisq(score.y^2,df = 1, low = F)

  # pValue unique
  pval <- apply(cbind(pvalu.x,pvalu.y), 1, max)^2
  # qValue by FDR
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F)$qval
  # confounders
  conf <- mod$U
  conf1 <- conf[(nrow(mod$U)-length(X) + 1) : nrow(mod$U),]

  # mediation
  if (mediation == T) {
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = cbind(conf1,conf.known), sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)
  }
  else {
    media <- NULL
  }

  return(list(pval  = pval,
              qval  = qval,
              conf  = conf,
              conf1 = conf1,
              score = cbind(score.x,score.y),
              mod   = mod,
              media = media))
}

#' wrapBACO : Mediation EWAS and independent mediation analysis (no coufounders estimated)
#'
#' @param X exposure
#' @param Y outcome
#' @param M methylation matrix
#' @param conf.known known confondeurs
#' @param mediation TRUE or FALSE, if TRUE : independant mediation analysis
#' @param GLMwrap result of wrapGLM
#' @param top If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()
#' @param sims Number of MC simulations in the mediation model
#'
#' @return pval, qval, score for each CpG
#' @return media : if mediation = T, result of wrapMed() function
#' @return mod : whole bacon model
#' 
#' @details 
#' 
#' The BACON model does not consider confusers but it aims to reduce the bias and inflation
#'  of a statistical test that in our case is a glm. 
#' Note that this method can be combined with a method estimating confondeurs.
#'
#' @examples
#' # Simulate data :
#' simu <- r_mediation(100, 500, 5)
#' # do a mediation EWAS and the independant mediation analysis
#' mod <- wrapBACO(simu$X, simu$Y, mediation = T)
#' # quick result of mediation EWAS
#' plot(-log10(mod$pval))
#' points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
#'
#' # or with a previous glm
#' mod1 <- wrapGLM(simu$X, simu$Y, mediation = F)
#' mod2 <- wrapBACO(simu$X, simu$Y, mediation = T, GLMwrap = mod1)
#' @export
wrapBACO <- function(X,Y,M,conf.known=NULL,mediation = T,top = 50,GLMwrap = NULL, sims = 100) {
  # GLM
  if (is.null(GLMwrap) == F){
    coeff <- GLMwrap
  }
  else {
    print("glm wait ...")
    coeff <- wrapGLM(X,Y,M,conf = conf.known)
  }


  # meta BACON model
  bc  <- bacon::bacon(effectsizes = coeff$es.bac,standarderrors = coeff$se.bac)
  mbc <- bacon::meta(bc)
  pval <- bacon::pval(mbc)
  pval <- pval[,1:2]
  pval <- apply(pval, 1, max)^2

  # qValue by FDR with meta pValue
  qval <- fdrtool::fdrtool(pval,statistic = "pvalue", verbose=F,plot=F)$qval
  # score unique (third column)
  score <- bacon::tstat(mbc)
  score <- score[,1:2]

  print("association done")

  # mediation
  if (mediation == T) {
    print("mediation wait ...")
    ord <- order(pval)[1:top]
    med <- M[,ord]

    ACME.eff <- NULL
    ACME.p.med <- NULL

    for (i in 1:top) {
      m <- wrapMed(X,Y,mediator = med[,i],CONF = conf.known, sims = sims)

      ACME.eff <- c(ACME.eff, m$ACME.eff)
      ACME.p.med <- c(ACME.p.med, m$ACME.p.med)
    }

    media <- data.frame(ACME.eff,ACME.p.med,cpg = ord)
  }
  else {
    media <- NULL
  }

  return(list(pval  = pval,
              qval  = qval,
              score = score,
              mod   = mbc,
              media = media))
}





