### Implementation of decorana in R

#' Decorana: R implementation
#'
#' Similar to Fortran implementation in vegan, but all in R.
#'
#' @export
`decorana` <-
    function(x, iweigh = 0, iresc = 0, ira = 0, mk = 26, short = 0,
             before = NULL, after = NULL, ...)
{
    .NotYetImplemented(iweigh)
    .NotYetImplemented(iresc)
    .NotYetImplemented(short)
    .NotYetImplemented(before)
    .NotYetImplemented(after)
    ## constants
    NAXES <- 4
    ## initialize: in vegan & standard CA style
    x <- vegan:::initCA(x)
    aidot <- attr(x, "RW")
    adotj <- attr(x, "CW")
    ## orthogonal CA: just do it and return
    if (ira == 1)
        return(orthoCA(x, aidot = aidot, adotj = adotj, naxes = NAXES))
}

#' Orthogonal Correspondence Analysis
#'
#' Orthogonal correspondence analysis is performed via svd of
#' CA-initialized data.
#'
#' @param x initCA-initialized data.
#' @param aidot,adotj Row and column weights summing up to 1.
#' @param naxes Number of axes.
#'
#' @return Orthogonal CA scaled like in Decorana.
#'
## not exported
`orthoCA` <-
    function(x, aidot, adotj, naxes)
{
    m <- svd(x, nu = naxes, nv = naxes)
    lambda <- m$d[seq_len(naxes)]^2
    rproj <- (m$u / sqrt(aidot)) %*% diag(sqrt(lambda/(1-lambda)), nrow = naxes)
    cproj <- (m$v / sqrt(adotj)) %*% diag(sqrt(1/(1-lambda)), nrow = naxes)
    list(evals = lambda, rproj = rproj, cproj = cproj)
}

#' Detrending
#'
#' Detrends axis x on mk segments of x1
#'
#' @param x Axis to be detrended.
#' @param aidot Weights of x.
#' @param x1 Axis along which x is detrended.
#' @param mk Number of segments on x1.
#'
#' @return Detrended values.
#'
## Not exported
`detrend` <-
    function(x, aidot, x1, mk)
{
    x1 <- cut(x1, mk)
    ## pad segments with zeros to buffer ends
    z <- c(0, 0, tapply(aidot*x, x1, sum), 0, 0)
    zn <- c(0, 0, tapply(aidot, x1, sum), 0, 0)
    z[!is.finite(z)] <- 0
    zn[!is.finite(zn)] <- 0
    ## mean by blocks of three segments
    z <- (z + c(0, z[-length(z)]) + c(z[-1],0))/
        pmax(c(zn + c(0, zn[-length(zn)]) + c(zn[-1], 0)), .Machine$double.eps)
    ## mean of blocks of three by starting position
    z <- (z + c(0, z[-length(z)]) + c(z[-1],0))/3
    x - z[as.numeric(x1)+2]
}
