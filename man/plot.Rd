\name{plot}

\alias{plot.JMbayes}

\title{MCMC Diagnostics for Joint Models}

\description{
  Produces MCMC diagnostics plots based on the \bold{coda} package.
}

\usage{
\method{plot}{JMbayes}(x, which = c("trace", "autocorr"), 
    param = c("betas", "tau", "inv.D", "gammas", "alphas", "Dalphas", 
        "sigma.t", "Bs.gammas"), ask = TRUE, \dots)
}

\arguments{
  \item{x}{an object inheriting from class \code{JMbayes}.}
  \item{which}{which types of plots to produce.}
  \item{param}{for which parameter to produce the MCMC diagnostic plots; default is for all parameters.}
  \item{ask}{logical, if \code{TRUE} the user is asked for input, before a new figure is drawn.}
  \item{\dots}{additional arguments; currently none is used.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.
}

\seealso{\code{\link{jointModelBayes}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)

# survival regression fit
fitSURV <- coxph(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# joint model fit, under the (default) Weibull model
fitJOINT <- jointModelBayes(fitLME, fitSURV, timeVar = "year")

plot(fitJOINT)
}
}

\keyword{methods}
