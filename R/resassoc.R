#' Residual Species Associations in Constrained Ordination
#'
#' Function finds residual species association after fitting the
#' constraints (and conditions) in constrained ordination
#' (\code{\link[vegan]{cca}}, \code{\link[vegan]{rda}}).
#'
#' @details
#'
#' The residual unconstrained axes can be used to assess if there is
#' any important unexplained variation in constrained ordination. This
#' approach was proposed already in the first papers on constrained
#' ordination (ter Braak 1986). However, this has been rarely done,
#' mainly because there are no easily available tools to assess the
#' unconstrained axes. Bayesian analysis of species communities
#' (Tikhonov et al. 2020) is conceptually similar to constrained
#' ordination. There the random effects are analysed \emph{via}
#' correlated species responses that are interpreted to be caused by
#' unknown environmental variables or residual inter-species
#' association. The random effects are found as Bayesian latent
#' factors which are conceptually similar to residual axes in
#' constrained ordination (although these are maximum likelihood
#' principal components). It is customary represent these residual
#' correlations as ordered correlation matrix (Tikhonov et al
#' 2017). This function provides similar tools for constrained
#' correpsondence analysis (CCA) and redundancy analysis (RDA) as
#' perfomed in \CRANpkg{vegan} functions \code{\link[vegan]{cca}} and
#' \code{\link[vegan]{rda}}.
#'
#' In Bayesian framework, we only try to extract a low number of
#' latent factors needed to describe correlated random variation. In
#' constrained ordination, all residual variation is accounted for by
#' residual axes. For meaningful analysis, we should only look at some
#' first axes which may reflect systematic unexplained variation. With
#' all axes, the residual correlations are usually nearly zero and
#' systematic variation is masked by non-systematic noise. Therefore
#' we should only have a look at some first axes. Users can either
#' specify the number of axes used, and the default is to use
#' brokenstick distribution (\CRANpkg{vegan} function
#' (\code{\link[vegan]{bstick}}) to assess the number of axes that
#' have surprisingly high eigenvalues and may reflect unaccounted
#' systematic variation. The function returns inter-species
#' associations as correlation-like numbers, where values nearly zero
#' mean non-associated species. These values can be plotted in
#' correlation plot or analysed directly. The correlations are
#' directly found from the residual ordination axes.
#'
#' It is an attractive idea to interpret inter-species association to
#' describe ecological interactions among species (Tikhonov et
#' el. 2017). This indeed is possible, but there are other alternative
#' explanations, such as missing environmental variables or modelling
#' errors (model formula, scaling of species, scaling of environmental
#' variables, assumptions on the shape of species response etc.). You
#' should be cautious in interpreting the results (and the same
#' applies also to the analogous Bayesian model).
#'
#'
#' @references
#'
#' ter Braak, C.J.F. (1986) Canonical correspondence analysis: a new
#' eigenvector technique for multivariate direct gradient
#' analysis. \emph{Ecology} \strong{67,} 1167--1179.
#'
#' Tikhonov, G., Abrego, N., Dunson, D. and Ovaskainen, O. (2017) Using
#' joint species distribution models for evaluating how
#' species-to-species associations depend on the environmental
#' context. \emph{Methods in Ecology and Evolution} \strong{8,} 443-
#' 452.
#'
#' Tikhonov et al. (2020) Joint species distribution modelling with
#' the R-package Hmsc. \emph{Methods in Ecology and Evolution}
#' \strong{11,} 442--447.
#'
#' @seealso \CRANpkg{Hmsc} has analogous function
#'     \code{computeAssociations}. However, in \pkg{Hmsc} the
#'     associations are based on Bayesian latent factors which are an
#'     essential natural component of the analysis whereas here we
#'     apply post-analysis tricks to extract something similar.
#'
#' @param x Constrained ordination result from
#'     \code{\link[vegan]{cca}} or \code{\link[vegan]{rda}}.
#' @param rank Number of unconstrained ordination axes that are used
#'     to compute the species associations. If missing, the rank is
#'     chosen by brokenstick distribution
#'     (\code{\link[vegan]{bstick}}).
#'
#' @return Correlation-like association matrix between species with
#'     added argument \code{"rank"} of the unconstrained data used to
#'     compute associations.
#'
#' @examples
#' library(vegan)
#' data(dune, dune.env)
#' mod <- cca(dune ~ A1 + Moisture, dune.env)
#' resassoc.cca(mod)
#'
#' @importFrom stats cov2cor
#' @importFrom vegan bstick
#'
#' @export
`resassoc.cca` <-
    function(x, rank)
{
    if (!inherits(x, "cca"))
        stop("input must be a constrained ordination object from vegan")
    if (inherits(x, c("dbrda","capscale")))
        stop("distance-based methods do not have species scores")
    if (is.null(x$CA$eig))
        stop("object 'x' does not have residual ordination")
    ## use bstick to get the rank if not given
    if (missing(rank)) {
        bs <- bstick(x$CA$rank, x$CA$tot.chi)
        rank <- min(which(bs > x$CA$eig)) - 1
    }
    ## Residual associations can be estimated only when rank > 1, else
    ## return identity matrix
    if (rank < 2) {
        r <- diag(NROW(x$CA$v))
        dimnames(r) <- list(rownames(x$CA$v), rownames(x$CA$v))
    } else {
        eig <- sqrt(x$CA$eig[seq_len(rank)])
        x <- x$CA$v[, seq_len(rank)] %*% diag(eig, nrow=rank)
        r <- cov2cor(tcrossprod(x))
    }
    attr(r, "rank") <- rank
    r
}
