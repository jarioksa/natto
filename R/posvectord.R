#' Position Vector and Species Vector Ordination
#'
#' Position Vector Ordination (PVO) finds ordination axes that go
#' through sample points in species space and maximize the projections
#' of other points onto them. These are similar to Principal
#' Components. Species Vector Ordination (SVO) finds a set of species
#' that maximize the projections of other species on species
#' axes. This also can be similar to Principal Components and can be
#' used to select a subset of species that cover most of the variance
#' in the data.
#'
#' Position Vector Ordination (PVO, function \code{posvectord}) is a
#' simple educational method that resembles Principal Component
#' Analysis, PCA (Orlóci 1966, 1973a). Function \code{posvectord}
#' positions vectors from the data centroid (origin) to sample points
#' in species space and evaluates the projections of other sample
#' vectors on these positioned axes, and selects the one with highest
#' total projection as the ordination axis. Then the effects of that
#' selected axis are removed from the covariance matrix, zeroing the
#' row and column of selected sampling unit, and the process is
#' repeated. Principal Components maximize this projection, but PVO
#' axes can be close to the Principal Component, in particular in the
#' large data sets with many observations. The method was proposed as
#' a computationally light approximation to PCA suitable for the
#' computers of 1960s (Orlóci 1966). Now it can only be regarded as
#' historically interesting, and also as an educational tool in
#' introducing PCA.
#'
#' Species Vector Ordination (SVO; function \code{spvectord}) is
#' similar, but it picks the species vector that maximizes the species
#' projected onto that vector (Orlóci 1973b).  The method picks the
#' species or variable that explains largest proportion of variance
#' and uses the centred values of this variable as the axis. Then it
#' orthogonalizes all species or variables to that selected axis and
#' repeats the selection process. PCA finds a linear combination of
#' species or variables that minimizes the residual variance, but in
#' SVO only one species or variable is used. The axes are named by the
#' species or variables, and the axis scores are the centred
#' (residual) observed values. The resulting plot has orthogonalized
#' set of species as axes. SVO can be similar to PCA, in particular
#' when some few species contribute most of the the total
#' variance. The method only has educational use in explaining PCA in
#' species space. Orlóci (1973b) suggesed SVO as a method of selecting
#' a subset of species or variables that contributes most of the
#' variance in the data. The current \code{spvectord} function can
#' also be used for this purpose, although it is mainly geared for
#' ordination. \CRANpkg{dave} package has function \code{orank}
#' specifically written for variable or species selection using the
#' same algorithm.
#'
#' @encoding UTF-8
#'
#' @param x Input data.
#' @param scale Scale variables to unit variance.
#' @param \dots Other arguments (passed to \code{\link[vegan]{ordiplot}}).
#'
#' @return \code{posvectord} returns an object of class
#'     \code{"posvectord"}, and \code{spvectord} returns an object of
#'     class \code{"spvectord"} that inherits from
#'     \code{"posvectord"}. Both result objects have the following
#'     elements:
#'
#' \itemize{
#'    \item \code{points}: The ordination scores. In SVO, these are
#'     named by the species (variable) the axes is based on, the
#'     numerical scores are the centred (residual) values of
#'     observed data. In PVO, they are called as PVO, and they are
#'     scaled similarly as in PCA and can reconstitute
#'     (a low-rank approximation of) covariances.
#'    \item \code{totvar}: The total variance in the input data.
#'    \item \code{eig}: Eigenvalues of axes.
#'}
#'
#' @seealso \code{\link{polarord}} (Polar Ordination) is a similar
#'     educational tool to approximate Principal Coordinates Analysis
#'     (PCoA).
#'
#' @references
#'
#' Orlóci, L. (1966) Geometric models in ecology. I. The theory and
#' application of some ordination methods. \emph{J. Ecol.} 54: 193--215.
#'
#' Orlóci, L. (1973a) Ordination by resemblance matrices. In:
#' R. H. Whittaker (ed.) \emph{Ordination and Classification of
#' Communities. Handbook of Vegetation Science} 5: 249--286.
#'
#' Orlóci, L. (1973b) Ranking characters by dispersion
#' criterion. \emph{Nature} 244: 371--373.
#'
#' @examples
#' data(spurn)
#' spvectord(spurn)
#' m <- posvectord(spurn)
#' m
#' plot(m)
#' if (require(vegan, quietly = TRUE)) {
#' plot(procrustes(rda(spurn), m, choices=1:2))
#' }

#' @rdname posvectord
#' @export
`posvectord` <-
    function(x, scale = FALSE)
{
    EPS <- sqrt(.Machine$double.eps)
    ## centre by species, crossproduct for sites
    x <- scale(x, scale = scale)
    ## variance if unscaled
    if (!scale)
        x <- x/sqrt(nrow(x)-1)
    u <- matrix(NA, nrow(x), ncol(x))
    eig <- numeric(ncol(x))
    names(eig) <- colnames(u) <- paste0("PVO", seq_len(ncol(x)))
    rownames(u) <- rownames(x)
    xx <- tcrossprod(x)
    totvar <- sum(diag(xx))
    for (i in 1:ncol(u)) {
        a <- xx/sqrt(pmax(diag(xx), EPS))
        crit <- rowSums(a^2)
        if (max(crit, na.rm = TRUE) < EPS)
            break
        take <- which.max(crit)
        u[,i] <- xx[take,]/sqrt(xx[take, take])
        eig[i] <- crit[take]
        xx <- xx - tcrossprod(u[,i])
        ## xx[take,] and xx[,take] should be zero, but may be, say, -2.2e-16
        xx[take,] <- 0
        xx[,take] <- 0
    }
    nit <- !is.na(colSums(u))
    out <- list("points" = u[, nit, drop=FALSE],
                "totvar" = totvar,
                "eig" = eig[nit])
    class(out) <- "posvectord"
    out
}
#'
#' @importFrom stats cov
#' @rdname posvectord
#' @export
`spvectord` <-
    function (x, scale = FALSE)
{
    x <- scale(x, scale = scale)
    u <- matrix(NA, nrow(x), ncol(x))
    eig <- numeric(ncol(x))
    colnames(u) <- paste0("Spec", seq_len(ncol(u)))
    rownames(u) <- rownames(x)
    xx <- cov(x)
    totvar <- sum(diag(xx))
    for(i in 1:ncol(u)) {
        ## off-diagonal elements of xx[i,j] are r[i,j]*s[i]*s[j],
        ## diagonal s^2[i], and we divide to give off-diagonal
        ## r[i,j]*s[j], diagonal s[i]
        xx <- xx/sqrt(diag(xx))
        ## explained variance r[i,j]^2 * s^2[j]
        crit <- rowSums(xx^2)
        if (max(crit, na.rm = TRUE) < sqrt(.Machine$double.eps))
            break
        take <- which.max(crit)
        eig[i] <- crit[take]
        u[,i] <- x[,take]
        colnames(u)[i] <- names(take)
        ## orthogonalize: x will be residuals
        x <- qr.resid(qr(cbind(1, u[,i])), x)
        ## x[,take] should now be exactly 0, but may be EPS (magnitude
        ## 1e-16), and then EPS/EPS in crit > 0
        x[,take] <- 0
        xx <- cov(x)
    }
    names(eig) <- colnames(u)
    nit <- !is.na(colSums(u))
    out <- list("points" = u[, nit, drop=FALSE],
                "totvar" = totvar,
                "eig" = eig[nit])
    class(out) <- c("spvectord", "posvectord")
    out
}

#' @export
`print.posvectord` <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    if (inherits(x, "spvectord"))
        cat("\nSpecies ")
    else
        cat("\nPosition ")
    cat("Vector Ordination\n\n")
    cat("Total Inertia ", round(x$totvar, digits), "\n\n")
    cat("Eigenvalues of Vectors:\n")
    print(zapsmall(x$eig), digits=digits)
    invisible(x)
}

#' @importFrom vegan scores ordiplot
#' @rdname posvectord
#' @export
`plot.posvectord` <-
    function(x, ...)
{
    ordiplot(x, display = "sites", ...)
}
