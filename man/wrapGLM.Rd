% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDMA.R
\name{wrapGLM}
\alias{wrapGLM}
\title{wrapGLM : naive glm function for mediation EWAS and independent mediation analysis (no estimated confondeurs)}
\usage{
wrapGLM(X, Y, M, conf.known = NULL, mediation = F, top = 20,
  family = "binomial", sims = 100)
}
\arguments{
\item{X}{exposure}

\item{Y}{outcome}

\item{M}{methylation matrix}

\item{conf.known}{known confondeurs}

\item{mediation}{TRUE or FALSE, if TRUE : independant mediation analysis}

\item{top}{If mediation is true, then "top" is the number of CpGs analyzed by wrapMed()}

\item{family}{"binomial" or "gaussian" : family function for GLM}

\item{sims}{Number of MC simulations in the mediation model}
}
\value{
pval, qval, score, standard errors and effect sizes for each CpG
}
\description{
wrapGLM : naive glm function for mediation EWAS and independent mediation analysis (no estimated confondeurs)
}
\details{
This function uses a generalized linear model "glm" to perform association tests of the mediation EWAS. 
Note that no confusers are estimated in this function.
}
\examples{
# Simulate data :
simu <- r_mediation(100, 500, 5)
# do a mediation EWAS and the independant mediation analysis
mod <- wrapGLM(simu$X, simu$Y, mediation = T)
# quick result of mediation EWAS
plot(-log10(mod$pval))
points(simu$mediators, -log10(mod$pval[simu$mediators]),col = 2)
}
