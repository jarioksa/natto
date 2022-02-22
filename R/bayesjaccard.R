#' Community Dissimilarity as Expected or Sampled Beta Variate
#'
#' Jaccard dissimilarity for presence/absence data can be seen as the
#' mode of Beta distributed variate. The function calculates the
#' dissimilarity as a random sample of Beta distribution, or
#' alternatively as its expected value or mode.
#'
#' @details
#'
#' In often-used 2x2 contingency table notation, \eqn{a} is the number
#' of species shared by two compared communities, and \eqn{b} and
#' \eqn{c} are the numbers of species occurring only in one of the
#' compared communities. Assuming uniform prior in (0, 1), the species
#' will \dQuote{sample} the dissimilarity of two communities with
#' posterior Beta(\eqn{b+c+1}, \eqn{a+1}). This will have mode
#' \eqn{(b+c)/(a+b+c)} and expected value \eqn{(b+c+1)/(a+b+c+2)}. The
#' mode is directly the Jaccard dissimilarity for binary
#' (presence/absence) data, and the expected value adds 1 in numerator
#' and 2 in denominator from the prior. The importance of prior will
#' diminish when the number of species \eqn{a+b+c} grows, but it will
#' protect from claiming complete identity or complete difference when
#' we only have a few species.
#'
#' Function \code{bayesjaccard} estimates Jaccard dissimilarity as
#' Beta-distributed random variate, and will return random sample from
#' its posterior distribution. It can also return the expected value
#' or the mode which are constant in function calls.
#'
#' The function is intended to be used to assess the effect of random
#' sampling variation in community analysis. The \pkg{natto} package
#' provided three examples of its application: \code{\link{clsupport}}
#' finds the \dQuote{support} of branches in hierarchic clustering by
#' function \code{\link{hclust}}, \code{\link{bjNMDS}} the variation
#' of ordination scores in non-metric multidimensional scaling by
#' functions \code{\link[vegan]{metaMDS}} and
#' \code{\link[vegan]{monoMDS}}, and \code{\link{bjdbrda}} the
#' variation of ordination scores of constrained component of
#' distance-based RDA by function \code{\link[vegan]{dbrda}}. All
#' these functions find the basic result from the expected value of
#' the dissimilarity, and produce a set of random samples from the
#' Beta distribution to assess the variation in the result.
#'
#' @seealso Function is a wrapper to \code{\link[vegan]{designdist}}.
#'
#' @return Dissimilarity object inheriting from classes \code{"dist"}
#'     and \code{"designdist"}.
#'
#' @examples
#'
#' data(spurn)
#' ## the effect of prior
#' plot(bayesjaccard(spurn, "mode"), bayesjaccard(spurn, "expected"))
#' abline(0, 1, col=2)
#' ## one random sample of dissimilarities
#' plot(bayesjaccard(spurn, "expected"), bayesjaccard(spurn), asp=1)
#' abline(0, 1, col=2)
#'
#' @param x Community data, will be treated as binary.
#' @param method Dissimilarity as a random sample from Beta
#'     distribution or as its expected value or mode.
#' @importFrom stats rbeta
#' @importFrom vegan designdist
#' @rdname bayesjaccard
#' @export
`bayesjaccard` <-
    function(x, method = c("rbeta", "expected", "mode"))
{
    method <- match.arg(method)
    switch(method,
           "expected" =
               designdist(x, "(b+c+1)/(a+b+c+2)", terms = "binary", abcd=TRUE),
           "mode" =
               designdist(x, "(b+c)/(a+b+c)", terms = "binary", abcd = TRUE),
           "rbeta" =
               designdist(x, "rbeta(length(a), b+c+1, a+1)", terms = "binary",
                          abcd = TRUE)
           )
}

#' Support of Branches of Hierarchical Clustering from Beta Jaccard
#'
#' Function performs hierarchic clustering (\code{\link{hclust}}) with
#' the expected value of Jaccard dissimilarity as Beta random variate,
#' and compares the stability of each cluster branch in random samples
#' from Beta distribution.
#'
#' @details
#'
#' Function is basically a graphical tool that plots an
#' \code{\link{hclust}} dendrogram and adds the count of similar
#' clusters in randomized cluster dendrograms, but it can also be
#' called only for the numeric result without plot.
#'
#' The count of similar clusters in trees is based on randomized form
#' Jaccard Beta dissimilarities, and is called here
#' \dQuote{support}. The support is calculated alternatively as the
#' number of exactly similar branches in randomized trees or as a sum
#' of maximum Jaccard similarities in observed and randomized trees
#' (\dQuote{softmatch}). The exact matches are very sensitive to
#' single wandering sampling units, and often give very low
#' \dQuote{support} for large classes, whereas Jaccard-based
#' \dQuote{support} may give too high values for the smallest
#' classes. For exact matches the \dQuote{support} is given as count,
#' and for soft maches (Jaccard-based) as a rounded integer per 1000.
#'
#' @seealso \code{\link{hclust}}, \code{\link{bayesjaccard}}.
#'
#' @return Usually called to draw a plot, but will return the
#'     \dQuote{support}; the values are returned in the order of
#'     merges in the agglomerative hierarchic clustering. Function
#'     \code{clsets} returns a list with indices of items (sampling
#'     units) of each branch in their merge order.
#'
#' @examples
#'
#' data(spurn)
#' clsupport(spurn) # exact match
#' clsupport(spurn, softmatch = TRUE)
#'
#'
#' @param x Community data, will be treated as binary.
#' @param n Number of random samples from Beta distribution.
#' @param method Clustering method, passed to \code{\link{hclust}}.
#' @param softmatch Estimate cluster similarity as Jaccard similarity
#'     or as a hard exact match between two clusters.
#' @param plot Plot the \code{\link{hclust}} dendrogram from expected
#'     Jaccard dissimilarity with each brand labelled with support.
#' @param ... Other parameters; passed to \code{\link{plot.hclust}}.
#'
#' @importFrom stats hclust rbeta
#' @importFrom vegan designdist ordilabel
#' @rdname clsupport
#' @export
`clsupport` <-
    function (x, n=1000, method="average", softmatch = FALSE, plot = TRUE, ...)
{
    h0 <- hclust(bayesjaccard(x, method="expected"), method)
    m0 <- clsets(h0)
    supp <- s <- numeric(nrow(h0$merge)-1)
    for(i in seq_len(n)) {
        d <- bayesjaccard(x)
        h <- hclust(d, method)
        m <- clsets(h)
        if (softmatch) {
            ## Find the Jaccard similarities of each m0 & m sets and
            ## take the highest as the value of "soft match". aabc
            ## (a+b + a+c) is the sum of species numbers of two sets,
            ## and abc the size of union of sets (number of units in
            ## pooled sets, a+b+c).
            aabc <- outer(sapply(m0, length), sapply(m, length), "+")
            abc <- matrix(0, length(m), length(m))
            for(i in seq_along(m))
                abc[,i] <- sapply(m0, function(x) length(unique(c(m[[i]], x))))
            abc <- aabc/abc - 1 # Jaccard similarity = a/abc
            s <- apply(abc, 1, max)
        } else {
        for (i in seq_along(supp))
            s[i] <- any(sapply(m, function(x) all(identical(m0[[i]], x))))
        }
        supp <- supp + s
    }
    if (plot) {
        ## with softmax supp is a float
        if (softmatch) {
            supp <- round(1000 * supp/n)
        }
        plot(h0, ...)
        ordilabel(h0, "internal", labels=supp)
    }
    supp
}

#' @param hclus \code{\link{hclust}} result object.
#' @rdname clsupport
#' @export
`clsets` <-
    function(hclus)
{
    m <- hclus$merge
    memb <- list()
    sets <- list()
    for (j in seq_len(nrow(m) - 1)) {
        for (k in 1:2)
            sets[[k]] <- if (m[j, k] < 0)
                             -m[j,k]
                         else
                             memb[[m[j,k]]]
        memb[[j]] <- sort(c(sets[[1]], sets[[2]]))
    }
    memb
}

#############
### NMDS  ###
#############

#' Multiple NMDS Ordinations from Beta Distributed Jaccard Dissimilarity
#'
#' Function performs \code{\link[vegan]{metaMDS}} on expected Beta
#' Jaccard dissimilarity with multiple starts and then a number of
#' \code{\link[vegan]{monoMDS}} runs on random Beta Jaccard
#' dissimilarities starting from the expected configuration, and
#' Procrustes rotates the random solutions to the expected one.
#'
#' @details
#'
#' Current function is designed for visual inspection of random
#' variation of ordination, and no numerical analysis functions are
#' (yet) available. The \code{plot} can show the scatter of points as
#' convex hulls or ellipsoid hulls containing a given proportion of
#' random points, or as stars that connect the expected point to
#' randomized ones. With option \code{"wedge"} the convex hull is
#' extended to include the observed point if necessary. See
#' \code{\link{bjpolygon}} and \code{\link{bjstars}} for technical
#' details.
#'
#' Missing functionality includes adding species scores, adding fitted
#' environmental vectors and factors and numeric summaries of the
#' randomization results.
#'
#' @return
#'
#' A \code{\link[vegan]{metaMDS}} result object amended with item
#' \code{rscores} that is a three-dimensional array of coordinates
#' from random \code{\link{bayesjaccard}} ordinations.
#'
#' @examples
#'
#' data(spurn)
#' m <- bjNMDS(spurn)
#' plot(m, keep = 0.607, col = "skyblue")

#' @param x Community data; will be treated as binary.
#' @param n Number of random samples of Beta Distribution.
#' @param trymax Maximum number of random starts in
#'     \code{\link[vegan]{metaMDS}}.
#' @param maxit,smin,sfgrmin,sratmax Convergence parameters in
#'     \code{\link[vegan]{monoMDS}}.
#' @param parallel Number of parallel tries in \code{\link[vegan]{metaMDS}}.
#' @param trace Trace iterations in \code{\link[vegan]{metaMDS}}.
#' @param ... Other parameters passed to functions.

#' @importFrom stats fitted
#' @importFrom vegan metaMDS monoMDS procrustes
#' @rdname bjNMDS
#' @export
`bjNMDS` <-
    function(x, n = 100, trymax = 500, maxit = 1000, smin = 1e-4,
             sfgrmin = 1e-7, sratmax = 0.999999, parallel = 2, trace=FALSE,
             ...)
{
    ## Expected ordination
    d0 <- bayesjaccard(x, method="expected")
    m0 <- metaMDS(d0, trymax = trymax, maxit = maxit, smin = smin,
                  sfgrmin = sfgrmin, sratmax = sratmax, parallel = parallel,
                  trace = trace, ...)
    ## random Jaccard ordinations
    rscore <- array(dim = c(dim(m0$points), n))
    for (i in seq_len(n)) {
        d <- bayesjaccard(x)
        m <- monoMDS(d, m0$points, maxit = maxit, smin = smin,
                     sfgrmin = sfgrmin, sratmax = sratmax, ...)
        rscore[,,i] <- fitted(procrustes(m0$points, m$points), truemean=FALSE)
    }
    m0$nsamp <- n
    m0$rscores <- rscore
    m0$call <- match.call()
    m0$distance <- "binary bayesjaccard"
    m0$data <- deparse(substitute(x))
    class(m0) <- c("bjnmds", class(m0))
    m0
}

#' @export
`print.bjnmds` <-
    function(x, ...)
{
    cat("NMDS based on BayesJaccard dissimilarity with", x$nsamp, "samples\n")
    NextMethod("print", x, ...)
    cat("Information refers to the expected ordination: samples will differ\n\n")
}

#' @param x Community data to be analysed and treated as binary
#'     (\code{bjNMDS}) or object to be plotted (\code{plot}).
#' @param choices Axes to be plotted.
#' @param kind Shape to be plotted to show the scatter of coordinates
#'     of sampled distances; see \code{\link{bjpolygon}} for details.
#' @param keep Proportion of points to be enclosed by shape; see
#'     \code{\link{peelhull}} and \code{\link{peelellipse}}.
#' @param type Type of the plot: \code{"t"}ext, \code{"p"}oints or
#'     \code{"n"}one.
#'
#' @rdname bjNMDS
#' @export
`plot.bjnmds` <-
    function(x, choices = 1:2, kind = c("hull", "ellipse", "wedge", "star"),
             keep = 0.9, type = "t", ...)
{
    kind <- match.arg(kind)
    x0 <- x$point[, choices, drop=FALSE]
    xarr <- x$rscores[, choices, , drop = FALSE]
    plot(x0, type = "n", asp = 1, xlab = paste0("NMDS", choices[1]),
         ylab = paste0("NMDS", choices[2]))
    switch(kind,
           "hull" = bjpolygon(xarr, x0, kind = "hull", observed = FALSE,
                              keep = keep, type = type, ...),
           "ellipse" = bjpolygon(xarr, x0, kind = "ellipse",
                                 keep = keep, type = type, ...),
           "wedge" = bjpolygon(xarr, x0, kind = "hull", observed = TRUE,
                               keep = keep, type = type, ...),
           "star" = bjstars(xarr, x0, keep = keep, type = type, ...)
           )
}

#' Shapes to Display Scatter of Points in Multiple Ordinations
#'
#' These are low-level functions that are used by \code{plot}
#' functions to draw shapes enclosing a specified proportion of
#' randomized points, or connecting certain proportion of them to the
#' reference points. The functions are not usually called directly by
#' users, but the plots can be modified with described arguments.
#'
#' @details
#'
#' Functions require output of three-dimensional array \code{xarr} of
#' randomized scores, where the last dimension is for the random
#' samples, and a two-dimensional matrix of references scores
#' \code{x0}.
#'
#' Function \code{bjpolygon} uses \code{\link{peelhull}} or
#' \code{\link{peelellipse}} to find a convex hull
#' (\code{\link{chull}}) or an ellipsoid hull
#' (\code{\link[cluster]{ellipsoidhull}}) containing \code{keep}
#' proportion of points. The reference coordinates can also be added
#' to the plot, either as point or a text label
#' (\code{\link{ordiellipse}}) and optionally connected to the shape
#' centre (not to the centre of points). With argument \code{observed}
#' the convex hull can be extended to include the reference point when
#' that is outside the hull; this is used to draw wedges for arrows
#' using the origin as the reference point.
#'
#' Function \code{bjstars} connects the reference point to \code{keep}
#' proportion of closest points. If the reference point is not
#' specified, the centre of points prior to selection will be used.
#'
#' @param xarr 3-D array of coordinates of sampling units times two
#'     axes by random samples.
#' @param x0 2-D array of sampling units times two axes treated as
#'     constant for all samples in \code{xarr}.
#' @param keep Proportion of points enclosed in the shape; passed to
#'     \code{\link{peelhull}} or \code{\link{peelellipse}}.
#' @param kind Shape is either a convex hull (\code{\link{peelhull}})
#'     or ellipse \code{\link{peelellipse}} enclosing \code{keep}
#'     proportion of \code{xarr} points.
#' @param criterion Criterion to remove extreme points from the hull
#'     in \code{\link{peelhull}} and \code{\link{peelellipse}}.
#' @param linetopoint Draw line from the centre of the shape to the
#'     coordinates in \code{x0}.
#' @param col Colour of the shape; can be a vector of colours.
#' @param alpha Transparency of shapes; 0 is completely transparent,
#'     and 1 is non-transparent.
#' @param observed After finding the convex hull, extend hull to
#'     enclose fixed point \code{x0}.
#' @param type Mark coordinate of \code{x0} using \code{"t"}ext,
#'     \code{"p"}oint or \code{"n"}one.
#' @param ... Other parameters passed to to marker of \code{type}.
#'
#' @return Functions return invisibly \code{NULL}.
#'
#' @seealso \code{\link{peelhull}}, \code{\link{peelellipse}},
#'     \code{\link{polygon}} for basic functions and
#'     \code{\link{plot.bjnmds}} and \code{\link{plot.bjnmds}} and
#'     \code{\link{plot.bjdbrda}} for programmatic interface.
#'
#' @importFrom graphics polygon
#' @importFrom grDevices adjustcolor
#' @rdname bjpolygon
#' @export
`bjpolygon` <-
    function(xarr, x0, keep = 0.9, kind = c("hull", "ellipse"),
             criterion = "area", linetopoint = TRUE, col="gray",
             alpha = 0.3, observed = TRUE, type = c("t", "p", "n"),
             ...)
{
    kind <- match.arg(kind)
    dims <- dim(xarr)
    nobs <- dims[1]
    nsam <- dims[3]
    ## handle colours
    if (alpha > 1)
        alpha <- alpha / 255
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length = nobs)
    ## draw polygons
    for (i in seq_len(nobs)) {
        poly <- switch(
            kind,
            "hull" = peelhull(t(xarr[i,,]), keep = keep,
                              criterion = criterion),
            "ellipse" = peelellipse(t(xarr[i,,]), keep = keep)
        )
        if (kind == "hull" && observed)
            poly <- peelhull(rbind(poly, x0[i,]), keep = 1)
        polygon(poly, col = adjustcolor(col[i], alpha.f = alpha), border = NA)
        if (linetopoint) {
            cnt <- attr(poly, "centre")
            ## line is non-transparent
            segments(x0[i,1], x0[i,2], cnt[1], cnt[2], col = col[i])
        }
    }
    switch(type,
           "n" = NULL,
           "p" = points(x0, col = col, ...),
           "t" = ordilabel(x0, ...)
           )
    invisible()
}

#' @importFrom graphics points segments
#' @importFrom vegan ordilabel
#' @rdname bjpolygon
#' @export
`bjstars` <-
    function(xarr, x0, keep = 0.9, col="gray", type = c("t", "p", "n"), ...)
{
    nobs <- nrow(x0)
    nsam <- dim(xarr)[3]
    nkept <- ceiling(nsam * keep)
    if (nkept == nsam)
        kept <- rep(TRUE, nsam)
    col <- rep(col, length = nobs)
    if (missing(x0))
        x0 <- colMeans(xarr)
    for (i in seq_len(nobs)) {
        if (nkept < nsam) {
            dist <- colSums((xarr[i,,] - x0[i,])^2)
            kept <- rank(dist) <= nkept
        }
        segments(x0[i,1], x0[i,2], xarr[i,1,kept], xarr[i,2,kept], col = col[i],  ...)
    }
    switch(type,
           "n" = invisible(),
           "p" = points(x0, ...),
           "t" = ordilabel(x0, ...)
           )
    invisible()
}


#############
### dbRDA ###
#############

### Somewhat trickier to implement than NMDS. (1) dbRDA is an
### eigenvector method, and if we rotate, eigenvalues would change,
### and we need to skip rotation, but fix the axis reflection.
### (2) The "expected" model can have (one) negative eigenvalue, but
### rbeta models can have several and variable numbers of negative
### eigenvalues. (3) Probably we should not allow any adjustment
### against negative eigenvalues as these are data-set dependent and
### destroy the beauty of the distance (expect sqrt.dis?).

#' Multiple dbRDA from Beta Distributed Jaccard Dissimilarity
#'
#' Function performs distance-based RDA on expected Beta Jaccard
#' dissimilarities, and then reruns the analysis on given number of
#' random Beta Jaccard dissimilarities.
#'
#' @details
#'
#' Function works only with constrained and partial constrained
#' ordination, and only saves the randomized results of constrained
#' ordination. To maintain the eigenvalues, function does not rotate
#' the random results to the reference ordination, but it will reflect
#' axes with reversed signs. The eigenvalues and magnitudes of
#' axiswise correlations to the reference ordination can be inspected
#' with \code{boxplot}.
#'
#' The display type in \code{plot} can be specified independently to
#' each type of score (or the score can be skipped). The default
#' display can be modified using a similarly named list of graphical
#' argument values that will replace the defaults. For instance, to
#' display linear combination scores as convex hull, use \code{lc =
#' "hull"}, and modify its parameters with list \code{lc.par}.
#'
#' @return Function returns \code{\link{dbrda}} result object amended
#'     with item \code{BayesJaccard} that is a list of randomized
#'     results with following items for each random sample:
#' \itemize{
#'   \item \code{tot.chi} total inertia
#'   \item \code{eig} eigenvalues of constrained axes
#'   \item \code{r} absolute value of correlation coefficient of
#'     randomized axis and the reference axis
#'   \item \code{u, wa, biplot, centroids} ordination scores for linear
#'     combination scores, weighted averages scores, biplot scores and centroid
#'     of factor constraints; these are unscaled and are usually accessed with
#'     \code{\link{scores.bjdbrda}} for appropriate scaling
#' }
#'
#' @examples
#' if (require(vegan)) {
#' data(dune, dune.env)
#' m <- bjdbrda(dune ~ A1 + Management + Moisture, dune.env)
#' ## eigenvalues and correlations showing axis stability
#' boxplot(m)
#' boxplot(m, kind = "correlation")
#' ## Default plot
#' plot(m)
#' ## modify plot
#' plot(m, wa = "ellipse",
#'    wa.par=list(col = dune.env$Management, keep = 0.607))
#' }
#'
#' @param formula,data Model definition of type \code{Y ~ Var1 + Var2,
#'     data = X}, where \code{Y} is dependent community data (handled
#'     as binary data), \code{Var1} and \code{Var2} are independent
#'     (explanatory) variables found in data frame \code{X}. See
#'     \code{\link[vegan]{dbrda}} for further information.
#' @param n Number of Jaccard dissimilarity matrices sampled from Beta
#'     distribution.
#' @param ... Other parameters passed to \code{\link[vegan]{dbrda}}.
#'
#' @importFrom vegan dbrda eigenvals
#' @rdname bjdbrda
#' @export
`bjdbrda` <-
    function(formula, data, n=100, ...)
{
    environment(formula) <- environment()
    m0 <- dbrda(formula = formula, data = data, distance="expected",
                dfun = bayesjaccard, ...)
    if (is.null(m0$CCA))
        stop("only implememented for constrained analysis")
    naxes <- length(m0$CCA$eig)
    ## sample, but first collect only eigenvalues
    tot.chi <- numeric(n)
    ev <- matrix(NA, nrow = n, ncol = naxes)
    r <- matrix(NA, nrow = n, ncol = naxes)
    u0 <- m0$CCA$u
    u <- array(dim = c(dim(u0), n))
    wa <- array(dim = c(dim(u0), n))
    bp <- array(dim = c(dim(m0$CCA$biplot), n))
    if (!is.null(m0$CCA$centroids))
        cn <- array(dim = c(dim(m0$CCA$centroids), n))
    else
        cn <- NULL
    for (i in 1:n) {
        m <- dbrda(formula, data, distance="rbeta", dfun = bayesjaccard, ...)
        tot.chi[i] <- m$tot.chi
        ev[i,] <- eigenvals(m, model = "constrained")
        ## Check reflected axes. 'r' is correlation because u,u0 have
        ## colSums 0 and colSums of squares 1
        nreal <- seq_len(ncol(m$CCA$u))
        r[i, nreal] <- colSums(u0[, nreal] * m$CCA$u)
        reflex <- diag(sign(r[i,nreal]))
        u[,nreal,i] <- m$CCA$u %*% reflex
        wa[,nreal,i] <- m$CCA$wa %*% reflex
        bp[,nreal,i] <- m$CCA$biplot %*% reflex
        ## u, wa & bp exist always, but centroids can be missing
        if (!is.null(cn))
            cn[,nreal,i] <- m$CCA$centroids %*% reflex
    }
    BJ <- list("tot.chi" = tot.chi, "eig" = ev, "r" = abs(r), "u" = u,
               "wa" = wa, "biplot" = bp, "centroids" = cn)
    m0$BayesJaccard <- BJ
    m0$call <- match.call()
    m0$inertia <- "squared binary bayesjaccard distance"
    class(m0) <- c("bjdbrda", class(m0))
    m0
}

#' @param x \code{bjdbrda} result object.
#' @param choices Selected ordination axes.
#' @param display,scaling,const Kind of scores, the scaling of scores
#'     (axes), and scaling constant with similar definitions as in
#'     \code{\link[vegan]{scores.rda}}.
#' @param expected Return scores of the expected ordination instead of
#'     ordination based on random samples from Jaccard dissimilarity.
#' @param ... Other parameters passed to functions; passed
#'     \code{\link[vegan]{dbrda}} in \code{bjdbrda}, ignored in
#'     \code{scores}.
#'
#' @importFrom vegan scores
#' @importFrom stats nobs
#' @rdname bjdbrda
#' @export
`scores.bjdbrda` <-
    function(x, choices = 1:2, display = c("wa", "lc", "bp", "cn"),
             scaling = "species", const,
             expected = TRUE, ...)
{
    if (!missing(const))
        .NotYetUsed(const)
    if (expected)
        return(NextMethod("scores", x, ...))
    display <- match.arg(display)
    nr <- nobs(x)
    n <- length(x$BayesJaccard$tot.chi)
    ## get scaling
    scales <- c("none", "sites", "species", "symmetric")
    if (!is.numeric(scaling)) {
        scaling <- match.arg(scaling, scales)
        scaling <- match(scaling, scales) - 1L
    }
    sco <- switch(display,
                  "wa" = x$BayesJaccard$wa[, choices, ],
                  "lc" = x$BayesJaccard$u[, choices, ],
                  "bp" = x$BayesJaccard$biplot[, choices, ],
                  "cn" = x$BayesJaccard$centroids[, choices, ]
                  )
    ## cycle through every sample for scaling
    for (i in seq_len(n)) {
        sumev <- x$BayesJaccard$tot.chi[i]
        const <- if(display == "bp") 1 else sqrt(sqrt((nr-1) * sumev))
        ev <- x$BayesJaccard$eig[i, choices]
        slambda <- switch(abs(scaling),
                          sqrt(ev/sumev),
                          rep(1, length(choices)),
                          sqrt(sqrt(ev/sumev))
                          )
        if (!is.null(slambda))
            sco[,,i] <- const * sweep(sco[,,i], 2, slambda, "*")
    }
    attr(sco, "score") <- display
    attr(sco, "scaling") <- scaling
    sco
}

#' @importFrom graphics points
#' @importFrom vegan ordilabel
`rdadraw` <-
    function(xarr, x0, kind, ...)
{
    kind <- match.arg(kind,
                      c("n", "p", "t", "hull", "ellipse", "wedge", "star"))
    switch(kind,
           "n" = NULL,
           "p" = points(x0, ...),
           "t" = ordilabel(x0, ...),
           "hull" = bjpolygon(xarr, x0, kind = "hull", observed = FALSE, ...),
           "ellipse" = bjpolygon(xarr, x0, kind = "ellipse", ...),
           "wedge" = bjpolygon(xarr, x0, kind = "hull", observed = TRUE,
                               ...),
           "star" = bjstars(xarr, x0, ...)
           )
}

#' @param wa,lc,cn,bp Display of corresponding scores. \code{"n"}
#'     skips the score, \code{"p"} and \code{"t"} use points or text
#'     for the expected score, and other shapes define
#'     \code{\link{bjpolygon}} or \code{\link{bjstars}} shape used for
#'     sampled random scores.
#' @param wa.par,lc.par,cn.par,bp.par List of arguments to modify the
#'     plotting parameters of the corresponding shape.
#' @param type Add \code{"t"}ext or \code{"p"}oint for the expected
#'     score to a shape of scatter of random shapes.
#'
#' @importFrom utils modifyList
#' @importFrom graphics arrows
#' @importFrom vegan ordiArrowMul ordiArrowTextXY ordilabel scores
#' @rdname bjdbrda
#' @export
`plot.bjdbrda` <-
    function(x, choices = 1:2, wa = "p", lc = "n", cn = "hull", bp = "wedge",
             wa.par = list(), lc.par = list(), cn.par = list(), bp.par = list(),
             scaling = "species", type = "t", ...)
{
    draw <- list(wa, lc, cn, bp) != "n"
    names(draw) <- c("wa","lc","cn","bp")
    display <- names(draw)[draw]
    ## NextMethod plot will give warnings of all extra parameters
    ## ("wa.par is not a graphical parameter" etc.) and therefore we
    ## turn warnings off
    op <- options(warn = -1)
    g <- NextMethod("plot", x, type = "n", display = display,
                    scaling = scaling)
    options(op) # back to user-defined old warn option
    ## Draw WA
    if (draw["wa"]) {
        xarr <- scores(x, choices = choices, display = "wa", scaling = scaling,
                       expected = FALSE)
        x0 <- scores(x, choices = choices, display = "wa", scaling = scaling,
                     expected = TRUE)
        def <- list(xarr = xarr, x0 = x0, kind = wa, col = "black", cex = 0.6)
        if (!wa %in% c("p", "t"))
            def <- modifyList(def,
                              list(col = "gray", alpha = 0.3, keep = 0.9,
                                   type = type))
        if (!is.null(wa.par))
            def <- modifyList(def, wa.par)
        do.call("rdadraw", def)
    }
    ## Draw LC
    if (draw["lc"]) {
        xarr <- scores(x, choices = choices, display = "lc", scaling = scaling,
                       expected = FALSE)
        x0 <- scores(x, choices = choices, display = "lc", scaling = scaling,
                     expected = TRUE)
        def <- list(xarr = xarr, x0 = x0, kind = lc, col = "darkgreen",
                    cex = 0.6)
        if (!lc %in% c("p", "t"))
            def <- modifyList(def,
                              list(alpha = 0.3, keep = 0.9, type = type))
        if (!is.null(lc.par))
            def <- modifyList(def, lc.par)
        do.call("rdadraw", def)
    }
    ## Draw centroids
    if (draw["cn"] && !is.null(g$centroids)) {
        xarr <- scores(x, choices = choices, display = "cn", scaling = scaling,
                       expected = FALSE)
        x0 <- scores(x, choices = choices, display = "cn", scaling = scaling,
                     expected = TRUE)
        def <- list(xarr = xarr, x0 = x0, kind = cn , col = "skyblue")
        if (!lc %in% c("p", "t"))
            def <- modifyList(def,
                              list(alpha = 0.3, keep = 0.9, type = type))
        if (!is.null(cn.par))
            def <- modifyList(def, cn.par)
        do.call("rdadraw", def)
    }
    ## Draw biplot arrows
    if (draw["bp"] && !is.null(g$biplot)) {
        arr <- ordiArrowMul(g$biplot)
        xarr <- arr * scores(x, choices = choices, display = "bp",
                             scaling = scaling, expected = FALSE)
        x0 <- arr * scores(x, choices = choices, display = "bp",
                           scaling = scaling, expected = TRUE)
        k <- rownames(x0) %in% rownames(g$biplot)
        xarr <- xarr[k,,, drop=FALSE]
        x0 <- x0[k,, drop=FALSE]
        orig <- matrix(0, nrow = nrow(x0), ncol=2)
        def <- list(xarr = xarr, x0 = x0, kind = bp, col = "blue")
        if (!lc %in% c("p", "t")) {
            def <- modifyList(def,
                              list(x0 = orig, alpha = 0.3, keep = 0.9,
                                   type = "n", lineto = FALSE))
            if (!is.null(bp.par))
                def <- modifyList(def, bp.par)
            do.call("rdadraw", def)
        }
        arrows(0, 0, x0[,1], x0[,2], col = def$col, length = 0.1)
        if (lc != "p")
            ordilabel(ordiArrowTextXY(x0, rescale=FALSE))
    }
}

## plot eigenvalues or axis correlations

#' @param kind Draw boxplots of eigenvalues or pairwise axis correlations
#'     for each random sample.
#' @param points,pch Add points of given color and shape to the
#'     eigenvalue boxplot.
#' @param xlab,ylab Change labelling of axes in boxplot.
#'
#' @importFrom graphics boxplot points
#' @rdname bjdbrda
#' @export
`boxplot.bjdbrda` <-
    function(x, kind = c("eigen", "correlation"), points = "red", pch = 16,
             xlab = "Axis", ylab, ...)
{
    kind <- match.arg(kind)
    if (kind == "eigen") {
        if (missing(ylab))
            ylab <- "Eigenvalue"
        bp <- boxplot(x$BayesJaccard$eig, xlab = xlab, ylab = ylab, ...)
        if (!is.na(points))
            points(x$CCA$eig, col = points, pch = pch, ...)
    } else if (kind == "correlation") {
        if (missing(ylab))
            ylab <- "Correlation"
        bp <- boxplot(x$BayesJaccard$r, xlab = xlab, ylab = ylab, ...)
        abline(h = 0, lty=3)
    }
    invisible(bp)
}
## print
##' @export
`print.bjdbrda` <-
    function(x, ...)
{
    cat("dbRDA based on BayesJaccard dissimilarity with",
        length(x$BayesJaccard$tot.chi), "samples\n\n")
    NextMethod("print", x, ...)
    cat("Information refers to the expected ordination: samples will differ\n\n")
}
