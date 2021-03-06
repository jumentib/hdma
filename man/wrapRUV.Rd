% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDMA.R
\name{wrapRUV}
\alias{wrapRUV}
\title{wrapRUV : Mediation EWAS and independent mediation analysis (coufounders estimated via RUV methods)}
\usage{
wrapRUV(X, Y, M, K, ctl = NULL, r.ctl = 0.2, meth = c("RUV2", "RUV4"),
  conf.known = NULL, mediation = T, top = 50, GLMwrap = NULL,
  sims = 100, ...)
}
\arguments{
\item{X}{exposure}

\item{Y}{outcome}

\item{M}{methylation matrix}

\item{K}{number of coufounders to estimate}

\item{ctl}{RUV methods need a control gene vector. if it is not known, it will be estimated by naive glms}

\item{r.ctl}{alpha risk for the estimation of the control gene vector}

\item{meth}{RUV combines two methods to estimate coufondeurs: "RUV2" and "RUV4"}

\item{conf.known}{known confondeurs}

\item{mediation}{TRUE or FALSE, if TRUE : independant mediation analysis}

\item{top}{If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()}

\item{GLMwrap}{result of wrapGLM to estimate the control gene}

\item{sims}{Number of MC simulations in the mediation model}

\item{...}{other options of the RUVs models}
}
\value{
pval, qval, score for each CpG

conf : coufounders estimated via RUV methods

media : if mediation = T, result of wrapMed() function

mod : whole RUV model

#' @details 

This function estimates the confusers via the RUV model. There are two sub-models of RUV to estimate the confondeurs, 
"RUV2" and "RUV4". Note that the RUV model uses control genes to estimate confounders. 
If you do not know it, it can be estimated via the wrapGLM function.
In addition, the RUV model use linear regression "lm" to perform associations tests of mediation EWAS.
}
\description{
wrapRUV : Mediation EWAS and independent mediation analysis (coufounders estimated via RUV methods)
}
\examples{
# Simulate data :
simu <- r_mediation(100, 500, 5)
# do a mediation EWAS and the independant mediation analysis
mod <- wrapRUV(simu$X, simu$Y, 5, mediation = T, meth ="RUV2)
# quick result of mediation EWAS
plot(-log10(mod$pval))
points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)

# or with a previous glm
mod1 <- wrapGLM(simu$X, simu$Y, mediation = F)
mod2 <- wrapRUV(simu$X, simu$Y, 5, mediation = T, meth ="RUV2, GLMwrap = mod1)
}
