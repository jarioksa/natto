### Proof-of-Concept function for constrained NMDS based on Bert van
### der Veen's idea as he explained to me in the IAVS meeting in Gijon
### in June 2026. Constrained NMDS is a monotonous (isometric)
### regression of observed dissimilarities on sample scores that are
### linear combinations of constraints (model matrix). This is a
### better idea than in natto::cNMDS which is a monotonous (isometric)
### regression of constrainted dissmilarities.

#' Constrained Non-metric Multidimensional Scaling
#'
#' Constrained non-metric Multidimensional scaling is a monotonous
#' (isometric) regression (code{\link[stats]{isoreg}}) of observed
#' (community) dissimilarities on dimensions that are linear
#' combinations of constraints.
#'
#' This is a proof-of-concept function based on the idea of Bert van
#' der Veen.
#'
#' @seealso \code(\link{cNMDS}} which provides an alternative and
#'     worse alternative that is bsed on monotonous (isomteric)
#'     regression of constrained dissimilarities
#'     \code{\link{distconstrain}}.
#'
#' @examples
#' data(mite, mite.env, package = "vegan")
#' mod <- caxNMDS(vegdist(mite) ~ WatrCont + SubsDens + Shrub + Topo,
#'   data = mite.env)
#' ordiplot(mod, display = "sites")
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
    u <- wcmdscale(D, k = k)
    B <- qr.coef(qr(mm), u)
    out <- optim(B, stress, D = D, mm = mm, k = k, method="L-BFGS-B")
    B <- matrix(out$par, ncol = k)
    rownames(B) <- colnames(mm)
    U <- mm %*% B
    ef <- envfit(U, df, permutations = 0)
    out <- list(formula = formula, coefficients = B, points = U, ef = ef)
    class(out) <- c("caxNMDS")
    out
}
