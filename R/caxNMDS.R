### Proof-of-Concept function for constrained NMDS based on Bert van
### der Veen's idea as he explained to me in the IAVS meeting in Gijon
### in June 2026. Constrained NMDS is a monotonous (isometric)
### regression of observed dissimilarities on sample scores that are
### linear combinations of constraints (model matrix). This is a
### better idea than in natto::cdisNMDS which is a monotonous (isometric)
### regression of constrainted dissmilarities.

#' Constrained Non-metric Multidimensional Scaling
#'
#' Constrained non-metric Multidimensional scaling is a monotonous
#' (isometric) regression (\code{\link[stats]{isoreg}}) of observed
#' (community) dissimilarities on dimensions that are linear
#' combinations of constraints.
#'
#' Ordinary non-constrained NMDS is a monotonous (isometric)
#' regression (\code{\link[stats]{isoreg}}) of observed (community)
#' dissimilarities on ordination scores. Constrained NMDS is similar,
#' but restricts the ordination scores to linear combinations of
#' constraining variables.
#'
#' Function uses standard \code{\link[stats]{optim}} to find the
#' regression coefficients of constraints. The regression coefficients
#' are the primary result of optimization, and ordination scores are
#' found applying regression coefficients on the model matrix of
#' constraints. The criterion variable to be minimized is
#' \code{stress} which is identical to stress in non-constrained NMDS
#' (\code{\link[vegan]{monoMDS}}, \code{\link[vegan]{metaMDS}} in
#' \pkg{vegan}). Function uses \dQuote{weak} ties, or allow breaking
#' tied values in monotonous regression.
#'
#' Function \code{caxNMDSengine} is the numerical core of the
#' method. It needs the input dissimilarities and model matrix
#' (without constant intercept). Function \code{caxNMDS} is more
#' user-friendly interface with model formula for constraints, and
#' also adds fitted vectors of continuous constraints and level
#' centroids of factor constraints. These are similar to biplot scores
#' and constraints in constrained ordination methods such as
#' \pkg{vegan} \code{\link[vegan]{rda}} and
#' \code{\link[vegan]{dbrda}}.
#'
#' Function \code{\link{cdisNMDS}} in this package provides an
#' alternative method based on constrained dissimilarities
#' (\code{\link{distconstrain}}) instead of constrained axes. The
#' current \code{caxNMDS} function seems to be more robust.
#'
#' This is a proof-of-concept function based on the idea of Bert van
#' der Veen.
#'
#' @seealso \code{\link{cdisNMDS}} which provides an alternative and
#'     worse alternative that is based on monotonous (isometric)
#'     regression of constrained dissimilarities
#'     \code{\link{distconstrain}}.
#'
#' @author Jari Oksanen and Bert van der Veen.
#'
#' @note Function \code{caxNMDSengine} can be used as \code{engine} in
#'     \code{\link[vegan]{metaMDS}} in \pkg{vegan} development version
#'     (2.8-0; not yet released).
#'
#' @return \code{caxNMDSengine} returns the \code{\link{optim}} result
#'     object and adds items \code{points} for ordination scores,
#'     \code{stress} for goodness of fit (\sQuote{stress}), and
#'     \code{coefficients} for a matrix of regression coefficients of
#'     constraints. Function \code{caxNMDS} returns these added items
#'     plus \code{formula}, \code{call} and fitted vectors and factor
#'     centroids of constraints in element \code{ef} with \pkg{vegan}
#'     function \code{\link[vegan]{envfit}}.
#'
#' @param formula Model formula where the left-hand-side is a distance
#'     structure for depenedent (community) dissimilarity and
#'     right-hand-side specifies the contraints.
#' @param data Data frame of constraints.
#' @param k Number of dimensions in NMDS.
#' @param method Optimization method used in \code{\link{optim}}.
#'
#' @examples
#' data(mite, mite.env, package = "vegan")
#' dis <- vegan::vegdist(mite)
#' mod <- caxNMDS(dis ~ WatrCont + SubsDens + Shrub + Topo,
#'   data = mite.env)
#' vegan::ordiplot(mod, display = "sites")
#' plot(mod$ef)
#' coef(mod)
#' mod$ef
#'

#' @importFrom vegan envfit
#' @export
`caxNMDS` <-
    function(formula, data, k = 2, method = "BFGS")
{
    ## Get data & response
    Trms <- delete.response(terms(formula, data = data))
    df <- model.frame(Trms, data = data)
    mm <- model.matrix(Trms, df)[,-1, drop=FALSE]
    mm <- scale(mm, scale=FALSE)
    D <- eval(formula[[2]], parent.frame(), environment(formula))
    ## optimize!
    sol <- caxNMDSengine(D = D, k = k, mm = mm, method = method)
    ## check & report optim result
    if (sol$convergence != 0)
        message("'optim' reported convergence issue ", sol$convergence,
                ": see ?optim")
    if (!is.null(sol$message))
        message(sol$message)
    ## output object
    ef <- envfit(sol$points, df, permutations = 0)
    out <- list(formula = formula, stress = sol$stress,
                coefficients = sol$coefficients, points = sol$points,
                ef = ef, call = match.call())
    class(out) <- "caxNMDS"
    out
}

#' @param D Dissimilarities.
#' @param u Initial configuration used for starting values in
#'     \code{\link{optim}}. Solution from metric scaling
#'     (\code{\link[vegan]{wcmdscale}}) is used if this is missing.
#' @param mm \code{\link{model.matrix}} of constraints.

#' @importFrom stats delete.response terms formula model.frame model.matrix
#' @importFrom stats dist isoreg optim
#' @importFrom vegan wcmdscale
#'
#' @rdname caxNMDS
#' @export
`caxNMDSengine` <-
    function(D, u, k, mm, method = "BFGS")
{
    ## stress
    stress <- function(B, D, mm, k) {
        B <- matrix(B, ncol = k)
        y <- dist(mm %*% B)
        ord <- order(D, y)
        s <- isoreg(y[ord])
        sqrt(sum((s$y - s$yf)^2)/sum(s$y^2))
    }
    ## analytic gradient of stress() wrt B (Kruskal 1964): the isotonic
    ## fit yf is already optimal for the current y, so its own
    ## sensitivity to a change in B drops out of the total derivative
    ## (envelope theorem) -- only the derivative of the Euclidean
    ## distance y wrt B remains
    stress_grad <- function(B, D, mm, k) {
        n <- nrow(mm)
        B <- matrix(B, ncol = k)
        Conf <- mm %*% B
        y <- dist(Conf)
        lt <- lower.tri(matrix(0, n, n))
        ord <- order(D, y)
        s <- isoreg(y[ord])
        yhat <- numeric(length(y))
        yhat[ord] <- s$yf

        Usum <- sum(y^2)
        Q <- sum((y - yhat)^2)/Usum
        dQ_dy <- 2*((y - yhat) - Q*y)/Usum

        ## spread the per-pair sensitivity over the points: each pair's
        ## sensitivity becomes a pull between the two points, and each
        ## point's total gradient is the sum of pulls from its partners
        W <- matrix(0, n, n)
        W[lt] <- dQ_dy/y
        W <- W + t(W)
        L <- diag(rowSums(W)) - W

        ## chain rule: stress = sqrt(Q), then Conf = mm %*% B
        as.vector(t(mm) %*% (L %*% Conf)) / (2*sqrt(Q))
    }
    if (missing(u) || is.null(u))
        u <- wcmdscale(D, k = k)
    B <- qr.coef(qr(mm), u)
    sol <- optim(B, stress, gr = stress_grad, D = D, mm = mm, k = k,
                 method = method)
    B <- matrix(sol$par, ncol = k)
    rownames(B) <- colnames(mm)
    U <- mm %*% B
    sol$points <- U
    sol$coefficients <- B
    sol$stress <- sol$value
    sol
}
