### Implementation of decorana in R

#' Decorana: R implementation
#'
#' Similar to Fortran implementation in vegan, but all in R.
#'
#' This function duplicates \CRANpkg{vegan} function
#' \code{\link[vegan]{decorana}}, and there is no need to use this
#' function for data analysis. The function serves two
#' purposes. Firstly, it is written in \R{} to allow easier inspection
#' of function than compiled code in C and Fortran in
#' \pkg{vegan}. Secondly, it is more hackable, and easier to develop
#' new features, change code or replace functionality than in the
#' compiled code. For instance, it would be trivial to add Detrended
#' Constrained Correspondence Analysis, but this would be impossible
#' without extensive changes in Fortran in the \pkg{vegan} function.
#'
#' Function is experimental, and at the moment it only implements
#' basic orthogonal and detrended correspondence analysis. It does not
#' implement many options of \code{\link[vegan]{decorana}}, most
#' importantly, rescaling is not yet implemented.
#'
#' @return Currently returns a list of elements \code{evals} of
#'     Decorana values, with \code{rproj} and \code{cproj} of scaled
#'     row and column scores.
#'
#' @examples
#' data(spurn)
#' str(mod <- rdecorana(spurn))
#' if (require(vegan)) {
#' ## compare to decorana
#' str(decorana(spurn, iresc = 0))
#' ## ordiplot works
#' ordiplot(mod, display="sites")
#' }
#'
#' @param x input data matrix.
#' @param iweigh Downweighting of rare species (0: no). Not yet
#'     implemented.
#' @param iresc Number of rescaling cycles (0: no rescaling). Not yet
#'     implemented.
#' @param ira Type of analyis (0: detrended, 1: orthogonal).
#' @param mk Number of segments in detrending.
#' @param short Shortest gradient to be rescaled. Not yet implemented.
#' @param before,after Definition of Hill's piecewise
#'     transformation. Not yet implemented.
#'
#' @export
`rdecorana` <-
    function(x, iweigh = 0, iresc = 0, ira = 0, mk = 26, short = 0,
             before = NULL, after = NULL)
{
    if (!missing(iweigh))
        .NotYetUsed("iweigh")
    if (!missing(iresc))
        .NotYetUsed("iresc")
    if (!missing(short))
        .NotYetUsed("short")
    if (!missing(before))
        .NotYetUsed("before")
    if (!missing(after))
        .NotYetUsed("after")
    ## constants
    NAXES <- 4
    EPS <- sqrt(.Machine$double.eps)
    CYCLES <- 1000
    ## initialize: in vegan & standard CA style
    x <- vegan:::initCA(x)
    aidot <- attr(x, "RW")
    adotj <- attr(x, "CW")
    ## orthogonal CA: just do it and return
    if (ira == 1)
        return(orthoCA(x, aidot = aidot, adotj = adotj, naxes = NAXES))
    ## ira == 0 and we go for detrended CA. First axis can be found
    ## directly via svd.
    rproj <- matrix(0, nrow(x), NAXES)
    cproj <- matrix(0, ncol(x), NAXES)
    evals <- numeric(NAXES)
    ## axis 1 is detrended
    sol <- svd(x, nu=1, nv=1)
    evals[1] <- sol$d[1]^2
    rproj[,1] <- sol$u / sqrt(aidot) * sqrt(evals[1]/(1-evals[1]))
    cproj[,1] <- sol$v / sqrt(adotj) * sqrt(1/(1-evals[1]))
    ## residual data
    x0 <- x
    x <- x  - tcrossprod(sol$u, sol$v) * sol$d[1]
    ## Go for the detrended axes
    for (axis in 2:NAXES) {
        ## svd of residual data
        sol <- svd(x, nu=1, nv=1)
        eig2 <- -1
        cycles <- 0
        ## Reciprocal averaging starting from eigenvector v
        repeat {
            sol <- transvu(sol$v[,1], rproj, x0, axis, aidot, adotj, mk)
            if (abs(eig2 - sol$d) < EPS)
                break
            eig2 <- sol$d
            if ((cycles <- cycles + 1) > CYCLES) {
                warning("no convergence on axis ", axis)
                break
            }
        }
        evals[axis] <- sol$d^2
        ## u is computed from v and includes eigenvalue
        u <- x0 %*% sol$v
        rproj[,axis] <- u / sqrt(aidot) * sqrt(1/(1-evals[axis]))
        cproj[,axis] <- sol$v / sqrt(adotj) * sqrt(1/(1-evals[axis]))
        ## residual matrix
        x <- x - tcrossprod(sol$u, sol$v)
    }
    list(evals = evals, rproj = rproj, cproj = cproj)
}

## Orthogonal Correspondence Analysis
##
## Orthogonal correspondence analysis is performed via svd of
## CA-initialized data.
##
## @param x initCA-initialized data.
## @param aidot,adotj Row and column weights summing up to 1.
## @param naxes Number of axes.
##
## @return Orthogonal CA scaled like in Decorana.
##
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

## transvu is modelled after trans subroutine in decorana.f. The
## decorana original is longer because it also implements
## orthogonalization which we do not need. Essentially the function
## only finds u from v, detrends u against previous axes, and finds
## new v. Their ratio is the eigenvalue of the step, and repeated call
## to this routine performs simple reciprocal averaging.  Function
## also returns normalized score vector v and eigenvalue-weighted u
## and the singular value (squareroot of the eigenvalue) so that the
## output has the same elements as svd().
#'
## not exported
`transvu` <-
    function(v, rproj, x, axis, aidot, adotj, mk)
{
    ## v should be centred and normalized: play safe
    cnt <- mean(sqrt(adotj) * v)
    v <- (v - cnt)
    v <- v / sqrt(sum(v^2))
    ## get u from v
    u <- x %*% v
    u <- u / sqrt(aidot)
    ## detrend
    if (axis > 1)
        for(k in c(seq_len(axis-1), rev(seq_len(axis-2))))
            u <- detrend(u, aidot, rproj[,k], mk)
    ## get back v
    v <- t(sqrt(aidot) * x) %*% u
    eig <- sqrt(sum(v^2))
    list(v = v / eig, u = sqrt(aidot) * u, d = sqrt(eig))
}

## Detrending
##
## Detrends axis x on mk segments of x1.
##
## @param x Axis to be detrended.
## @param aidot Weights of x.
## @param x1 Axis along which x is detrended.
## @param mk Number of segments on x1.
##
## @return Detrended values.
##
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
