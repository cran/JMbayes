\name{jointModelBayes}
\alias{jointModelBayes}

\title{
Joint Models for Longitudinal and Time-to-Event Data
}

\description{
Fits shared parameter joint models for longitudinal and survival outcomes under a Bayesian approach.
}

\usage{
jointModelBayes(lmeObject, survObject, timeVar, 
    survMod = c("weibull-PH", "spline-PH"), 
    param = c("td-value", "td-extra", "td-both", "shared-RE"), 
    robust = FALSE, robust.b = FALSE, df = 4, df.b = 4, lag = 0, init = NULL, 
    extraForm = NULL, priors = NULL, control = list(), ...)
}

\arguments{
  \item{lmeObject}{an object of class 'lme' fitted by function \code{lme()} from package \bold{nlme}.}
  \item{survObject}{an object of class 'coxph' fitted by function \code{coxph()} from package \bold{survival}.}
  \item{timeVar}{a character string indicating the time variable in the linear mixed effects model.}
  \item{survMod}{a character string specifying the type of the baseline risk function. See \bold{Details}.}
  \item{param}{a character string specifying the type of association structure between the longitudinal and 
        survival processes. See \bold{Details}.}
  \item{robust}{logical, if \code{TRUE} the linear mixed model is assumed to have a Student's-\eqn{t} errors distribution. See \bold{Details}.}
  \item{robust.b}{logical, if \code{TRUE} the random effects are assumed to have a Student's-\eqn{t} distribution. See \bold{Details}.}
  \item{df}{a numeric scalar denoting the number of degrees of freedom for the Student's-\eqn{t} error distribution for the mixed model.}
  \item{df.b}{a numeric scalar denoting the number of degrees of freedom for the Student's-\eqn{t} random-effects distribution.}
  \item{lag}{a numeric scalar denoting a lag effect in the time-dependent covariate represented by the mixed model; default is 0.}
  \item{init}{a named list of user-specified initial values:
    \describe{
        \item{betas}{the vector of fixed effects for the linear mixed effects model.}
        \item{tau}{the precision parameter from the linear mixed effects model (i.e., \eqn{\tau = 1/\sigma^2} with \eqn{\sigma} 
            denoting the error terms standard deviation).}
        \item{inv.D}{the inverse variance-covariance matrix of the random effects.}
        \item{b}{a matrix of random effects values.}
        \item{gammas}{the vector of baseline covariates for the survival model.}
        \item{alpha}{the association parameter(s).}
        \item{Dalpha}{the association parameter for the true slopes formulation.}
        \item{Bs.gammas}{the vector of spline coefficients for the spline-approximated baseline risk function; 
            specified only when \code{survMod = "spline-PH"}.}
        \item{sigma.t}{the shape parameter for the Weibull baseline risk function; specified only when \code{survMod = "weibull-PH"}.}
    }
    When this list of initial values does not contain some of these components or contains components 
    not of the appropriate length, then the default initial values are used instead.
  }
  \item{extraForm}{a list with components \code{fixed} a formula representing the fixed-effects part of the user-defined term, 
    \code{indFixed} a numeric vector indicating which fixed effects of \code{lmeObject} are involved in the user-defined term, 
    \code{random} a formula representing the random-effects part of the user-defined term, and \code{indRamdom} a numeric vector 
    indicating which random effects of \code{lmeObject} are involved in the user-defined term. Required only when 
    \code{param = "td-extra"} or \code{param = "td-both"}. See \bold{Examples}.}
  \item{priors}{a named list of user-specified prior parameters:
    \describe{
        \item{priorMean.betas}{the prior mean vector of the normal prior for the fixed effects of the linear mixed effects model.}
        \item{priorTau.betas}{the prior precision matrix of the normal prior for the fixed effects of the linear mixed effects model.}
        \item{priorA.tau}{the prior shape parameter of the Gamma prior for the precision parameter of the linear mixed effects model.}
        \item{priorB.tau}{the prior rate parameter of the Gamma prior for the precision parameter of the linear mixed effects model.}
        \item{priorMean.gammas}{the prior mean vector of the normal prior for the regression coefficients of the survival model.}
        \item{priorTau.gammas}{the prior precision matrix of the normal prior for the regression coefficients of the survival model.}
        \item{priorMean.alphas}{the prior mean vector of the normal prior for the association parameter in the survival model.}
        \item{priorTau.alphas}{the prior precision matrix of the normal prior for the association parameter in the survival model.}
        \item{priorMean.Dalphas}{the prior mean vector of the normal prior for the slope association parameter in the survival model.}
        \item{priorTau.Dalphas}{the prior precision matrix of the normal prior for the slope association parameter in the survival model.}
        \item{priorMean.Bs.gammas}{the prior mean vector of the normal prior for the spline coefficients of the baseline risk function;
            used only when \code{survMod = "spline-PH"}.}
        \item{priorTau.Bs.gammas}{the prior precision matrix of the normal prior for the spline coefficients of the baseline risk function; 
            used only when \code{survMod = "spline-PH"}.}
        \item{priorA.sigma.t}{the prior shape parameter of the Gamma prior for the shape parameter of the Weibull baseline risk function;
            used only when \code{survMod = "weibull-PH"}.}
        \item{priorB.sigma.t}{the prior rate parameter of the Gamma prior for the shape parameter of the Weibull baseline risk function;
            used only when \code{survMod = "weibull-PH"}.}
        \item{priorR.D}{the prior precision matrix of the Wishart prior for the precision matrix of the random effects.}
        \item{priorK.D}{the degrees of freedom of the Wishart prior for the precision matrix of the random effects.}
    }
  }
  \item{control}{a list of control values with components:
    \describe{        
        \item{program}{a character string specifying which program to use to fit the joint model. Available options are
            \code{WinBUGS}, \code{OpenBUGS} and \code{JAGS} (default).}
        \item{n.chains}{integer specifying the number of chains; default is 1.}
        \item{n.iter}{integer specifying the total number of iterations; default is 10000.}
        \item{n.burnin}{integer specifying how many of \code{n.iter} to discard as burn-in; default is 5000.}
        \item{n.thin}{integer specifying the thinning in each chain; default is 1.}
        \item{n.adapt}{integer specifying the number of iterations to use for adaptation; used only when \code{program = "JAGS"}
            and default is 1000.}
        \item{K}{integer used as the prior precision in the normal prior for the fixed effects, the regression coefficients of the survival submodel, 
            the association parameters, the slope association parameters, and in the spline coefficients; default is 100.}
        \item{C}{integer specifying the constant used in the zeros-trick; default is 5000.}
        \item{working.directory}{a character string specifying the working directory; default is the output of \code{getwd()}.}
        \item{bugs.directory}{a character string specifying the path to the \code{WinBUGS} executable; default is \code{"C:/Program Files/WinBUGS14/"}.}
        \item{openbugs.directory}{a character string specifying the path to the \code{OpenBUGS} executable; default is \code{NULL}.}
        \item{clearWD}{logical indicating whether to clean-up the working directory; default is \code{TRUE}.}
        \item{over.relax}{logical indicating whether over-relaxation should be employed; default is \code{TRUE}.}
        \item{bugs.seed}{integer specifying the seed for \code{WinBUGS} and \code{OpenBUGS}; default is 1.}
        \item{quiet}{a logical indicating whether to print information during the updating from \code{JAGS}; default is \code{FALSE}.}
        \item{knots}{a numeric vector of knots positions for the spline approximation of the log baseline risk function; default is \code{NULL}, which
            means that the knots are calculated based on the percentiles of the observed event times.}
        \item{lng.in.kn}{a numeric scalar indicating the number of knots to use (based on the percentiles); default is 5.}
        \item{ordSpline}{a numeric scalar setting the order of the spline function. This is the number of coefficients in each piecewise polynomial segment, 
	    thus a cubic spline has order 4; default is 4.}
        \item{ObsTimes.knots}{logical; if \code{TRUE} (default), the knots are set using the percentiles of the observed event times
            (i.e., including both true events and censored observations). If \code{FALSE}, the knots are set based on the percentiles of
            the true event times alone.}
    }
  }
  \item{\dots}{options passed to the \code{control} argument.}
}

\details{
Function \code{jointModelBayes} fits joint models for longitudinal and survival data under a Bayesian approach using WinBUGS, OpenBUGS or JAGS. 
For the longitudinal responses the linear mixed effects model represented by the \code{lmeObject} is assumed. If argument \code{robust} is set to 
\code{TRUE}, then linear mixed effects model fitted by \code{jointModelBayes} has the same fixed- and random-effects structure as in \code{lmeObject},
but the error terms distribution is assumed to be Student's-\eqn{t} with \code{df} degrees of freedom. If argument \code{robust.b} is set to 
\code{TRUE}, then the random-effects distribution is assumed to be a multivariate Student's-\eqn{t} with \code{df.b} degrees of freedom (when \code{FALSE} it is 
assumed to be a multivariate normal). For the survival times, let \eqn{w_i} denote the vector of baseline covariates in \code{survObject}, with associated parameter 
vector \eqn{\gamma}, \eqn{m_i(t)} the value of the longitudinal outcome at time point \eqn{t} as approximated by the linear mixed model (i.e., \eqn{m_i(t)} equals the 
fixed-effects part \code{+} random-effects part of the linear mixed effects model for sample unit \eqn{i}), \eqn{m_i'(t)} denotes an extra user-defined term to be
included in the linear predictor, and \eqn{\alpha} and \eqn{\alpha_d} vector of association parameters for \eqn{m_i(t)} and \eqn{m_i'(t)}, respectively. Then, \code{jointModelBayes}
assumes a relative risk model of the form \deqn{h_i(t) = h_0(t) \exp\{\gamma^\top w_i + f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d)\},} where the baseline risk
function is assumed to be of the Weibull form (i.e., \eqn{h_0(t) = \sigma_t t^{\sigma_t - 1}}) when \code{survMod = "weibull-PH"}, whereas for \code{survMod = "spline-PH"} is
approximated using splines, i.e., \deqn{\log h_0(t) = \sum_k \tilde\gamma_k B(t; \lambda),} with \eqn{B(.)} denoting a B-spline basis function, 
\eqn{\lambda} a vector of knots, and \eqn{\tilde \gamma_k} the corresponding splines coefficients.

Function \eqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d)} specifies the association structure between the two processes. In particular, for
\code{param = "td-value"}, \deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = \alpha m_i(t),} for \code{param = "td-extra"}, 
\deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = \alpha_d m_i'(t),} for \code{param = "td-both"}, 
\deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = \alpha m_i(t) + \alpha m_i(t),} and for \code{param = "shared-RE"}, 
\deqn{f(m_i(t), m_i'(t), b_i; \alpha, \alpha_d) = \alpha^\top b_i,} where \eqn{b_i} denotes the vector of random effects for 
the \eqn{i}th subject.
}

\value{
See \code{\link{JMbayesObject}} for the components of the fit.
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2000) Joint modelling of longitudinal measurements and event time data.
\emph{Biostatistics} \bold{1}, 465--480.

Hsieh, F., Tseng, Y.-K. and Wang, J.-L. (2006) Joint modeling of survival and longitudinal data: Likelihood 
approach revisited. \emph{Biometrics} \bold{62}, 1037--1043.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: With 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011) Dynamic predictions and prospective accuracy in joint models for longitudinal 
and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.

Tsiatis, A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. 
\emph{Statistica Sinica} \bold{14}, 809--834.

Wulfsohn, M. and Tsiatis, A. (1997) A joint model for survival and longitudinal data measured with error. 
\emph{Biometrics} \bold{53}, 330--339.
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}


\seealso{
\code{\link{coef.JMbayes}},
\code{\link{ranef.JMbayes}}, 
\code{\link{logLik.JMbayes}}, 
\code{\link{survfitJM}},
\code{\link{aucJM}},
\code{\link{dynCJM}},
\code{\link{prederrJM}}, 
\code{\link{predict.JMbayes}}
}

\examples{

\dontrun{
# A joint model for the AIDS dataset:
# First we fit the linear mixed model for the longitudinal measurements of
# sqrt CD4 cell counts
lmeFit.aids <- lme(CD4 ~ obstime * drug, random = ~ obstime | patient, data = aids)
# next we fit the Cox model for the time to death (including the 'x = TRUE' argument)
survFit.aids <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# the corresponding joint model is fitted by (the default is to assume 
# a Weibull relative risk, and the current value parameterization, and use
# JAGS to fit the model)
jointFit.aids <- jointModelBayes(lmeFit.aids, survFit.aids, timeVar = "obstime")
summary(jointFit.aids)

# A joint model for the PBC dataset:
# We first fit the linear mixed and Cox models. In the first we include 
# splines to model flexibly the subject-specific longitudinal trajectories
lmeFit.pbc <- lme(log(serBilir) ~ ns(year, 3),
    random = list(id = pdDiag(form = ~ ns(year, 3))), data = pbc2)
survFit.pbc <- coxph(Surv(years, status2) ~ 1, data = pbc2.id, x = TRUE)

# the corresponding joint model with a spline-approximated baseline risk
# function is fitted, using WinBUGS (you may have to adjust the 
# 'bugs.directory' control argument), by:
jointFit.pbc <- jointModelBayes(lmeFit.pbc, survFit.pbc, timeVar = "year", 
    survMod = "spline-PH", program = "winbugs", bugs.directory = "C:/WinBUGS/WinBUGS14/")
summary(jointFit.pbc)

# we update the joint model fitted for the PBC dataset by including
# the time-dependent slopes term. To achieve this we need to define 
# the 'extraForm' argument, in which we use function dns() to numerically
# compute the derivative of the natural cubic spline. In addition, we increase
# the number of MCMC iterations to 15000
dform = list(fixed = ~ 0 + dns(year, 3), random = ~ 0 + dns(year, 3),
    indFixed = 2:4, indRandom = 2:4)
jointFit.pbc2 <- update(jointFit.pbc, param = "td-both", extraForm = dform,
    n.iter = 15000)
summary(jointFit.pbc2)

# we fit the same model with the shared random effects formulation
jointFit.pbc3 <- update(jointFit.pbc, param = "shared-RE")
summary(jointFit.pbc3)
}
}

\keyword{multivariate}
\keyword{regression}
