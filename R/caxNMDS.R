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
#' This is a proof-of-concept function based on the idea of Bert van
#' der Veen.
#'
#' @seealso \code{\link{cdisNMDS}} which provides an alternative and
#'     worse alternative that is based on monotonous (isometric)
#'     regression of constrained dissimilarities
#'     \code{\link{distconstrain}}.
#'
#' @param formula Model formula where the left-hand-side is a distance
#'     structure for depenedent (community) dissimilarity and
#'     right-hand-side specifies the contraints.
#' @param data Data frame of constraints.
#' @param k Number of dimensions in NMDS.
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

#' @importFrom stats delete.response terms formula model.frame model.matrix
#' @importFrom stats dist isoreg optim
#' @importFrom vegan wcmdscale envfit
#' @export
`caxNMDS` <-
    function(formula, data, k = 2)
{
    ## Get data & response
    Trms <- delete.response(terms(formula, data = data))
    df <- model.frame(Trms, data = data)
    mm <- model.matrix(Trms, df)[,-1, drop=FALSE]
    mm <- scale(mm, scale=FALSE)
    D <- eval(formula[[2]], parent.frame(), environment(formula))
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
    u <- wcmdscale(D, k = k)
    B <- qr.coef(qr(mm), u)
    sol <- optim(B, stress, gr = stress_grad, D = D, mm = mm, k = k,
                 method="BFGS")
    ## check & report optim result
    if (sol$convergence != 0)
        message("'optim' reported convergence issue ", sol$convergence,
                ": see ?optim")
    if (!is.null(sol$message))
        message(sol$message)
    B <- matrix(sol$par, ncol = k)
    rownames(B) <- colnames(mm)
    U <- mm %*% B
    ef <- envfit(U, df, permutations = 0)
    out <- list(formula = formula, stress = sol$value, coefficients = B,
                points = U, ef = ef)
    class(out) <- c("caxNMDS")
    out
}
