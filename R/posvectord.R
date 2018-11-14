#' Position Vector Ordination
#'
#' Position Vector Ordination (PVO) is a simple educational method of
#' that resembles Principal Component Analysis (PCA). The method picks
#' the species or variable that explains largest proportion of
#' variance and uses the centred values of this variable as the
#' axis. Then it orthogonalizes all species or variables to that
#' selected axis and repeats the selection process. In PCA a linear
#' combination of species or variables is found to minimize the
#' residual variance, but in PVO only one species or variable is
#' used. The axes are named by the species or variables, and the axis
#' scores are the centred (residual) observed values.
#'
#' @encoding UTF-8
#'
#' @param x Input data.
#' @param scale Scale variables to unit variance.
#'
#' @return The function returns an object of class \code{"posvector"}
#'  with following elements:
#'  \itemize{
#'     \item \code{points}: The ordination scores. These are named by
#'       the species (variable) the axes is based on, the numerical
#'       scores are the centred (residual) values of observed data.
#'    \item \code{totvar}: The total variance in the input data.
#'    \item \code{eig}: Eigenvalues of axes.
#'  }
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
#' @rdname posvectord
#' @export
`posvectord` <-
    function(x, scale = FALSE)
{
    ## centre by species, crossproduct for sites
    x <- scale(x, scale = scale)
    ## variance if unscaled
    if (!scale)
        x <- x/sqrt(nrow(x)-1)
    u <- matrix(NA, nrow(x), ncol(x))
    eig <- numeric(ncol(x))
    names(eig) <- colnames(u) <- paste0("PVO", seq_len(ncol(x)))
    xx <- tcrossprod(x)
    totvar <- sum(diag(xx))
    for (i in 1:ncol(u)) {
        a <- xx/sqrt(pmax(diag(xx), .Machine$double.eps))
        crit <- rowSums(a^2)
        if (max(crit, na.rm = TRUE) < sqrt(.Machine$double.eps))
            break
        take <- which.max(crit)
        u[,i] <- xx[take,]/sqrt(xx[take, take])
        eig[i] <- crit[take]
        xx <- xx - tcrossprod(u[,i])
    }
    nit <- !is.na(colSums(u))
    out <- list("points" = u[, nit, drop=FALSE],
                "totvar" = totvar,
                "eig" = eig[nit])
    class(out) <- "posvector"
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
        x <- apply(x, 2, function(z) ortho(u[,i, drop=FALSE], z))
        xx <- cov(x)
    }
    names(eig) <- colnames(u)
    nit <- !is.na(colSums(u))
    out <- list("points" = u[, nit, drop=FALSE],
                "totvar" = totvar,
                "eig" = eig[nit])
    class(out) <- "posvector"
    out
}

