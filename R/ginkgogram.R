### This function is built upon core R stats:::plot.hclust function
### with the following information:

#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2016 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

#' Plot Dendrogram with Leaves as Fans
#'
#' Function is similar to standard \code{\link{hclust}}, but the
#' leaves are drawn as fans with base proportional to weights (or
#' sizes) of leaves.
#'
#' @param x An object of the type produced by \code{\link{hclust}}.
#' @param labels A character vector of labels for te leaves of the
#'     tree.
#' @param check Logical indicating if the \code{x} should be checked
#'     for validity.
#' @param axes,frame.plot,ann Logical flags as in
#'     \code{\link{plot.default}}.
#' @param main,sub,xlab,ylab Character strings to replace default
#'     annotation.
#' @param w Weights of leaves.
#' @param col Colour of the fill of leaves.
#' @param leaf Maximum height of leaf triangle. Leaf height is
#'     smaller of this value and the height of the branch, and the
#'     default value (\code{Inf}) always draws leaves to the
#'     root. Special value \code{TRUE} will take the smallest branch
#'     height as the maximum height of each leaf.
#' @param \dots Further graphical arguments
#'
#' @details Function \code{ginkgogram} has two new arguments to
#'     \code{\link{hclust}}: \code{w} for weights that give the base
#'     width of leaves, and \code{col} that gives the colour of the
#'     fill of leaves. The weights \code{w} must be numeric and
#'     non-negative, and zero weight guarantees that the leaf is drawn
#'     as a simple line similarly as in \code{\link{hclust}}.
#'
#' @examples
#' ## need vegan data sets
#' if (require(vegan)) {
#' data(dune, dune.phylodis, dune.taxon)
#' cl <- hclust(dune.phylodis)
#' ginkgogram(cl, w = colMeans(dune), col = factor(dune.taxon$Superorder))
#' ## limit heights of the polygons to the Cretaceous
#' ginkgogram(cl, w = colMeans(dune), col = factor(dune.taxon$Superorder),
#'     leaf = 65)
#' }

#' @importFrom grDevices dev.flush dev.hold
#' @importFrom graphics axis box plot.new polygon title
#' @export
`ginkgogram` <-
    function (x, labels = NULL, check = TRUE, axes = TRUE,
              frame.plot = FALSE, ann = TRUE, main = "Cluster Ginkgogram",
              sub = NULL, xlab = NULL, ylab = "Height", w, col,
              leaf = Inf, ...)
{
    hang <- -1
    ## check weights
    if (!missing(w) && !isTRUE(all(w >= 0)))
        stop("weights 'w' must be numeric and non-negative")
    merge <- x$merge
    if (check && !isTRUE(msg <- stats:::.validity.hclust(x, merge)))
        stop(msg)
    storage.mode(merge) <- "integer"
    n1 <- nrow(merge)
    n <- n1 + 1L
    if (missing(col))
        col <- 1
    col <- rep(col, length.out = n)
    ## get leafwidths and branch spreads
    if (missing(w))
        w <- rep(1, n)
    fanw <- w/sum(w)
    sprd <- pmax(1/n, fanw)
    span <- sum(sprd)
    fanw <- n1 * fanw / span
    fanw <- fanw - min(fanw)
    sprd <- n1 * sprd / span
    ordr <- cumsum(sprd[x$order]) + 1
    ordr <- ordr - sprd[x$order]/2
    #ordr <- ordr[order(x$order)]
    height <- as.double(x$height)
    labels <- if (missing(labels) || is.null(labels)) {
        as.character(if (is.null(x$labels)) seq_len(n) else x$labels)
    }
    else {
        if (is.logical(labels) && !labels)
            character(n)
        else as.character(labels)
    }
    dev.hold()
    on.exit(dev.flush())
    plot.new()
    graphics:::plotHclust(n1, merge, height, ordr[order(x$order)],
                          hang, labels, ...)
    ## draw polygons

    ## first we need polygon heights: either (1) to the root, or (2)
    ## up to a max certain value, or (3) as a special case, the lowest
    ## root (first height)
    if (isTRUE(leaf))
        leaf <- height[1] # the first is the lowest
    leaf <- pmin(height, leaf)
    for (i in seq_len(n1))
        for (j in 1:2) {
            if(merge[i,j] < 0) {
                k <- which(x$order == -merge[i,j])
                a <- fanw[x$order][k] / 2
                if (a > 0)
                    polygon(c(ordr[k]-a, ordr[k]+a, ordr[k]), c(0,0,leaf[i]),
                            col = col[-merge[i,j]])
            }
        }
    if (axes)
        axis(2, at = pretty(range(c(0, height))), ...)
    if (frame.plot)
        box(...)
    if (ann) {
        if (!is.null(cl <- x$call) && is.null(sub))
            sub <- paste0(deparse(cl[[1L]]), " (*, \"", x$method,
                "\")")
        if (is.null(xlab) && !is.null(cl))
            xlab <- deparse(cl[[2L]])
        title(main = main, sub = sub, xlab = xlab, ylab = ylab,
            ...)
    }
    invisible()
}
