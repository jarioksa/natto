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
    if (!missing(before))
        .NotYetUsed("before")
    if (!missing(after))
        .NotYetUsed("after")
    ## constants
    NAXES <- 4
    EPS <- sqrt(.Machine$double.eps)
    CYCLES <- 1000
    ## initialize: in vegan & standard CA style
    xorig <- x/sum(x)
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

## Original comments of Decorana detrnd:
##    starts with a vector x and detrends with respect to groups defined
##    by ix.  detrending is in blocks of 3 units at a time, and the
##    result calculated is the average of the 3 possible results that
##    can be obtained, corresponding to 3 possible starting positions
##    for the blocks of 3.

## NB: smoothing is done twice, and for equal weights zn effectively
## performs c(1,2,3,2,1)/3 smoothing, and defines smoothing window of
## 5 blocks.

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
#' @importFrom stats filter
##
## Not exported
`detrend` <-
    function(x, aidot, x1, mk)
{
    x1 <- cut(x1, mk)
    ## pad segments with zeros to buffer ends
    z <- c(0, 0, tapply(aidot*x, x1, sum, default = 0), 0, 0)
    zn <- c(0, 0, tapply(aidot, x1, sum, default = 0), 0, 0)
    ## mean z with weights zn by blocks of three
    z <- filter(z, c(1,1,1), sides=2)
    zn <- filter(zn, c(1,1,1), sides=2)
    z <- z / pmax(zn, .Machine$double.eps)
    ## mean of blocks of three
    z <- filter(z, c(1,1,1)/3, sides=2)
    x - z[as.numeric(x1)+2]
}

## The original comments of Decorana smooth:
##    takes a vector z and does (1,2,1)-smoothing until no blanks left
##    and then 2 more iterations of (1,2,1)-smoothing.  if no blanks to
##    begin with, then does 3 smoothings, i.e. effectively (1,6,15,20,
##    15,6,1)-smoothing.
##
## @param z vector to be smoothed.
##
#' @importFrom stats filter
##
## not exported
`smooth` <-
    function(z)
{
    kernel <- c(0.25, 0.5, 0.25) # (1,2,1)-smoothing
    mk <- length(z)
    idx <- seq_len(mk) + 1L
    repeat{
        z <- filter(c(z[1], z, z[mk]), kernel, sides=2)[idx]
        if (all(z > 0))
            break
    }
    z <- filter(c(z[1], z, z[mk]), kernel, sides=2)[idx]
    filter(c(z[1], z, z[mk]), kernel, sides=2)[idx]
}


## segment is the key function in rescaling and estimates the
## dispersion of species scores on each segment of the gradient (but
## does not do the segmenting). Original comments:
##   given an ordination (u, v), calculates numbers and summed
##   mean-square deviations in mk segments.  zn(k) is the number of
##   samples in segment k; zv(k) is the summed mean-square deviation.
##   (we aim to make zv, zn as nearly equal as possible.)

`segment` <-
    function(xorig, rproj, cproj, mk, aidot)
{
    ## collect statistics: these needs weights of original community
    ## data as species abundances are used as weights (and missing
    ## species do not contribute to site statistics).
    sqcorr <- rowSums(xorig^2)
    sumsq <- rowSums(xorig * outer(rproj, cproj, "-")^2)
    sqcorr <- pmin(sqcorr/aidot^2, 0.9999) # 0.9999 as in decorana.f
    sumsq <- sumsq/aidot
    axbit <- cut(rproj, mk)
    zv <- tapply(sumsq, axbit, sum, default = 0)
    zn <- tapply(1-sqcorr, axbit, sum, default = 0)
    list(zv = zv, zn = zn)
}

`stretch` <-
    function(x, xorig, rproj, cproj, aidot, short = 0)
{
    mk <- 20  # overrules user setting of mk
    z <- segment(xorig, rproj, cproj, mk, aidot)
    zv <- smooth(z$zv)
    zn <- smooth(z$zn)
    ## set within-sample square deviation to be 1
    sd <- sqrt(sum(zv/zn)/mk)
    rproj <- rproj/sd
    cproj <- cproj/sd
    along <- diff(range(rproj)) # length of axis
    if (along < short)
        return(list(rproj = rproj, cproj = cproj))
    ## new mk: not user-settable
    mk <- floor(5 * along) + 1L
    z <- segment(xorig, rproj, cproj, mk = mk, aidot)
    zv <- smooth(z$zv)
    zn <- smooth(z$zn)
    ## segment lenghts
    zv <- 1 / sqrt(0.2/along + zv/zn)
    zv <- zv * along/sum(zv)
    zn <- c(0, cumsum(zv))
    axbit <- along/mk
    ## species scores v are rescaled!
    iv <- trunc(cproj/axbit) + 1L
    iv <- pmin(pmax(iv, 1), mk)
    cproj <- zn[iv] + zv[iv] * (cproj/axbit - iv + 1)
    rproj <- (x %*% cproj) / aidot
    ## second pass
    mk <- 20
    z <- segment(xorig, rproj, cproj, mk, aidot)
    zv <- smooth(z$zv)
    zn <- smooth(z$zn)
    ## set within-sample square deviation to be 1
    sd <- sqrt(sum(zv/zn)/mk)
    rproj <- rproj/sd
    cproj <- cproj/sd
    list(rproj = rproj, cproj = cproj)
}

## Scale raw results from rescaling function stretch. In decorana this
## is part of eigy, but we want to pack that in a separate function to
## keep the main rdecorana() cleaner.

`postscale` <-
    function(z, xorig, aidot)
{
    axlong <- sqrt(sum(aidot * z$rproj^2))
    rproj <- z$rproj/axlong
    cproj <- z$cproj/axlong
    sumsq <- sum(xorig * outer(rproj, cproj, "-")^2)
    sd <- sqrt(sumsq)
    cproj <- cproj/sd
    list(rproj = rproj, cproj = cproj)
}
