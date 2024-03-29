% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gnlmod.R
\name{gnlmod}
\alias{gnlmod}
\alias{predict.gnlmod}
\alias{summary.gnlmod}
\title{Generalized Non-Linear Regression Analysis}
\usage{
gnlmod(formula, data, family = gaussian, wts, ...)

\method{predict}{gnlmod}(object, newdata, type = "response", ...)

\method{summary}{gnlmod}(object, dispersion = NULL, correlation = FALSE, ...)
}
\arguments{
\item{formula}{Model formula for a \code{\link{selfStart}}
\code{\link{nls}} model.}

\item{data}{Data for formula.}

\item{family}{Error distribution of the exponential
\code{\link{family}}.}

\item{wts}{Prior weights}

\item{\dots}{Other parameters passed to \code{\link{nlm}} that
performs the actual regression.}

\item{object}{\code{gnlmod} result object.}

\item{newdata}{New independent data for prediction.}

\item{type}{Type of predicted values; only \code{"response"} implemented.}

\item{dispersion}{The dispersion parameter for the family
used. Either a single numerical value of \code{NULL} (the
default), when it is inferred from the \code{object}.}

\item{correlation}{logical; if \code{TRUE}, the correlation matrix
of the estimated parameters is returned and printed.}
}
\value{
Function returns the \code{\link{nlm}} result object, but
    it is amended with \code{\link{glm}} object related to
    \code{\link{family}}.
}
\description{
Non-linear regression (\code{\link{nls}}) \code{\link{selfStart}}
model generalized to use exponential \code{\link{family}} error
distributions.
}
\details{
Function combines three \R{} functions for fitting non-linear
models with exponential family errors:
\enumerate{
  \item Non-linear regression is defined by
     \code{\link{selfStart}} like in \code{\link{nls}}. In addition to
     defining the nonlinear function, these also provide initial values and
     derivatives of model parameters.
  \item \code{\link{family}} functions are used to define the error
     distribution allowing use of other than least squares models. The
     \code{family} functions provide the log-likelihood function for fitting,
     and allows finding the partial derivatives of model parameters for the
     log-likelihood function (see McCullagh & Nelder 1989, p. 41).
  \item The log-likelihood function (with partial derivatives) is optimized
     with \code{\link{nlm}}.
}
The result is mostly similar to \code{\link{glm}} object, and can be
handled with many \code{glm} method functions, except those that assume
linear fit.
}
\examples{
## Compare to vegan species-area models
require(vegan) || stop("needs vegan")
data(sipoo, sipoo.map)
S <- specnumber(sipoo)
## Arrhenius model in vegan with least squares
nls(S ~ SSarrhenius(area, k, z), sipoo.map)
## Assume Poisson error
(mod <- gnlmod(S ~ SSarrhenius(area, k, z), sipoo.map, family=poisson))
## Arrhenius should be identical to Poisson glm with log-link
glm(S ~ log(area), sipoo.map, family=poisson) # (Intercept) = log(k)
## some method functions
fitted(mod)
predict(mod)
predict(mod, newdata = list(area = seq(1, 250, len=31)))
summary(mod)

}
\references{
McCullagh, P & Nelder, J.A. (1989) \emph{Generalized Linear Models}, 2nd ed.
   Chapman & Hall.
}
