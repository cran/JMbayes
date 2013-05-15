\name{JMbayesObject}

\alias{JMbayesObject}

\title{Fitted JMbayes Object}

\description{
An object returned by the \code{jointModelBayes} function, inheriting from class \code{JMbayes} and representing a fitted
joint model for longitudinal and time-to-event data. Objects of this class have methods for the generic functions 
\code{coef}, \code{confint}, \code{fixed.effects}, \code{logLik}, \code{plot}, \code{print}, 
\code{random.effects}, \code{summary}, and \code{vcov}.
}

\value{
  The following components must be included in a legitimate \code{JMbayes} object.
  \item{codaFit}{a list of class \code{mcmc.list} with the MCMC samples for each parameter, except from the random effects.}
  \item{postMeans}{a numeric matrix with the posterior means for the random effects of each subject.}
  \item{postVars}{a list with the posterior variance-covariance matrix for the random effects of each subject.}
  \item{modes}{a numeric vector with the approximate posterior modes for each parameter 
	(calculated via \code{density(..., bw = "nrd", adjust = 3, n = 1000))}.}
  \item{DIC}{the DIC value.}
  \item{pD}{the pD value.}
  \item{coefficients}{a list with the posterior means for the model's parameters. The components of this list are:
    \describe{
        \item{betas}{the vector of fixed effects for the linear mixed effects model.}
        \item{sigma}{the measurement error standard deviation for the linear mixed effects model.}
        \item{gammas}{the vector of baseline covariates for the survival model.}
        \item{alpha}{the association parameter(s).}
        \item{Dalpha}{the association parameter(s) corresponding to the slope of the true trajectory; returned only when
            \code{param = "td-extra"} or \code{param = "td-both"}.}
        \item{sigma.t}{the scale parameter for the Weibull survival model; returned only when
            \code{survMod = "weibull-PH"}.}
        \item{Bs.gammas}{the coefficients of the B-splines used to approximate the baseline hazard;
            returned only when \code{survMod = "spline-PH"}.}
        \item{D}{the variance-covariance matrix of the random effects.}
    }
  }
  \item{StErr}{a list with the same components as \code{coefficients} above with the standard deviations of the model's parameters.}
  \item{CIs}{a list with the same components as \code{coefficients} above with the 95\% confidence intervals of the model's parameters.}
  \item{vcov}{the variance-covariance matrix of the model's parameters based.}  
  \item{x}{a list with the design matrices for the longitudinal and event processes.}
  \item{y}{a list with the response vectors for the longitudinal and event processes.}
  \item{id}{the grouping vector for the longitudinal responses.}
  \item{times}{a numeric vector containing the time points at which the longitudinal responses were collected.}
  \item{data.id}{a \code{data.frame} containing the variables for the linear mixed effects model at the time of the
    event.}
  \item{survMod}{the value of the \code{survMod} argument.}
  \item{termsYx}{the \code{terms} component for the fixed-effects part of the \code{lmeObject}.}
  \item{termsYz}{the \code{terms} component for the random-effects part of the \code{lmeObject}.}
  \item{termsYx.deriv}{the \code{terms} component for the fixed-effects part of the slope term specified in argument \code{extraForm}.}
  \item{termsYz.deriv}{the \code{terms} component for the random-effects part of the slope term specified in argument \code{extraForm}.}
  \item{termsT}{the \code{terms} component of the \code{survObject}.}
  \item{formYx}{the formula for the fixed effects part of the longitudinal model.}
  \item{formYz}{the formula for the random effects part of the longitudinal model.}
  \item{formT}{the formula for the survival model.}
  \item{timeVar}{the value of the \code{timeVar} argument}
  \item{control}{the value of the \code{control} argument.}
  \item{param}{the value of the \code{param} argument.}
  \item{extraForm}{the value of the \code{extraForm} argument.}
  \item{robust}{the value of the \code{robust} argument.}
  \item{df}{the value of the \code{df} argument.}
  \item{priors}{a list with the specification of the prior distributions for the model's parameters. This has the same components as
    the \code{priors} argument of the \code{\link{jointModelBayes}} function.}
  \item{call}{the matched call.}
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{
    \code{\link{jointModelBayes}}
}

\keyword{multivariate}
\keyword{regression}
