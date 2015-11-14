#' Constrained and Residual Dissimilarities
#'
#' Function constrains dissimilarities by external variables, or
#' alternatively removes effects of constraining variables and returns
#' residual dissimilarities. The analysis is based on McArdle &
#' Anderson (2001), and the analysis of constrained dissimilarities is
#' equal to distance-based Redundancy Analysis
#' (\code{\link[vegan]{dbrda}}).
#'
#' @details Function uses the method of McArdle & Anderson (2001) to
#' constrain dissimilarities by external variables, or alternatively,
#' to find residual dissimilarities after constraints. With Euclidean
#' distances, the method is equal to performing linear regressions on
#' each column in the raw data and then calculationg the distances,
#' but works directly on distances. With other methods, there is no
#' similar direct connection to the raw data, but it is possible to
#' work with non-Euclidean metrics. The same basic method is used
#' within db-RDA (\code{\link[vegan]{dbrda}} in \CRANpkg{vegan}), but
#' this function exposes the internal calculations to users.
#'
#' Non-Euclidean indices can produce negative eigenvalues in
#' db-RDA. Would negative eigenvalues be produced, this function can
#' return negative squared distances resulting in \code{NaN} when
#' taking the square root. Db-RDA works with the internal presentation
#' of the dissimilarities, and its analysis does not suffer from the
#' imaginary distances, but these can ruin the analysis of
#' dissimilarities returned from this function.
#' 
#' @references McArdle, B.H. & Anderson, M.J. (2001). Fitting
#'   multivariate models to community data: a comment on distance-based
#'   redundancy analysis. \emph{Ecology} 82, 290--297.
#'
#' @importFrom stats delete.response terms model.frame model.matrix
#' formula as.dist
#'
#' @param formula The left-hand-side must be dissimilarities and the
#' right-hand-side should list the constraining variables.
#'
#' @param data Data frame containing the constrainging variables in
#' the \code{formula}.
#'
#' @param add an additive constant is added to the non-diagonal
#'    dissimilarities such that all \eqn{n-1} eigenvalues are
#'    non-negative. Alternatives are \code{"lingoes"} (default, also
#'    used with \code{TRUE}) and \code{"cailliez"} (which is the only
#'    alternative in \code{\link{cmdscale}}).
#' 
#' @param residuals Return residuals after constraints.
#'
#' @param squared Return squared dissimilarities instead of
#' dissimilarities. This allows handling negative squared distances by
#' the user instead of setting them \code{NaN}.
#' 
#' @export
`distconstrain` <-
    function(formula, data, add = FALSE, residuals = FALSE,
             squared = FALSE)
{
    ## evaluate data and get the model matrix
    if (missing(data))
        data <- .GlobalEnv
    Trms <- delete.response(terms(formula, data = data))
    df <- model.frame(Trms, data = data)
    mm <- model.matrix(Trms, df)[,-1, drop=FALSE]
    mm <- scale(mm, scale = FALSE) # centre
    ## extract dissimilarities (or the left-hand-side of formula) and
    ## perform the Gower standardization
    dis <- eval(formula[[2]])
    dis <- as.matrix(dis)
    ## handle add constant to make dis Euclidean
    if (is.logical(add) && isTRUE(add))
        add <- "lingoes"
    if (is.character(add)) {
        nd <- row(dis) != col(dis)
        add <- match.arg(add, c("lingoes", "cailliez"))
        if (add == "lingoes") {
            ac <- vegan:::addLingoes(dis)
            dis[nd] <- sqrt(dis[nd]^2 + 2 * ac)
        } else if (add == "cailliez") {
            ac <- vegan:::addCailliez(dis)
            dis[nd] <- dis[nd] + ac
        }
    } else {
        ac <- NA
    }
    dis <- -vegan:::GowerDblcen(dis^2)/2
    ## QR decomposition
    Q <- qr(mm)
    ## get fitted or residual dissimilarities
    if (residuals)
        dis <- qr.resid(Q, t(qr.resid(Q, dis)))
    else
        dis <- qr.fitted(Q, t(qr.fitted(Q, dis)))
    ## back to dissimilarities
    de <- diag(dis)
    dis <- -2 * dis + outer(de, de, "+")
    ## check and return
    dis <- as.dist(dis)
    if (!is.na(ac) && ac > 0) {
        attr(dis, "add") <- add
        attr(dis, "ac") <- ac
    }
    if (squared && any(dis < 0))
        warning("some dissimilarities were negative, use 'squared = TRUE'?")
    if (!squared)
        dis <- sqrt(dis)
    dis
}
