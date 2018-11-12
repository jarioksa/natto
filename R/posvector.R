#' Position Vector Ordination
#'
#' Position Vector Ordination (PVO) is a simple educational method of
#' that resembles Principal Component Analysis (PCA). The method pick
#' the species or variable that explains largest proportion of
#' variance and uses the centred values of this variable as the
#' axis. Then it orthogonalizes all species or variables to that
#' selected axis and repeats the selection process. In PCA a linear
#' combination of species or variables is found to minimize the
#' residual variance, but in PVO only one species or variable is
#' used. The axes are named by the species or variables, and the axis
#' scores are the centred (residual) observed values.
#'
#' @param x Input data.
#' @param scale Scale variables to unit variance.
#'
#' @importFrom stats cov
#'
#' @export
`posvector` <-
    function (x, scale = FALSE)
{
    x <- scale(x, scale = scale)
    u <- matrix(NA, nrow(x), ncol(x))
    colnames(u) <- paste0("Sp", seq_len(ncol(u)))
    for(i in 1:ncol(u)) {
        ## off-diagonal elements are r[i,j]*s[i]*s[j], diagonal s^2[i]
        xx <- cov(x)
        ## now off-diagonal are r[i,j]*s[j], diagonal s[i]
        xx <- xx/sqrt(diag(xx))
        ## explained variance r[i,j]^2 * s^[j]
        crit <- rowSums(xx^2)
        if (max(crit, na.rm = TRUE) < sqrt(.Machine$double.eps))
            break
        take <- which.max(crit)
        print(names(take))
        u[,i] <- x[,take]
        colnames(u)[i] <- names(take)
        ## orthogonalize: x will be residuals
        x <- apply(x, 2, function(z) ortho(u[,i, drop=FALSE], z))
    }
    u[, !is.na(colSums(u))]
}
