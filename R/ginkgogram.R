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

plot.hclust <-
function (x, labels = NULL, hang = 0.1, check = TRUE, axes = TRUE, 
    frame.plot = FALSE, ann = TRUE, main = "Cluster Dendrogram", 
    sub = NULL, xlab = NULL, ylab = "Height", ...) 
{
    merge <- x$merge
    if (check && !isTRUE(msg <- .validity.hclust(x, merge))) 
        stop(msg)
    storage.mode(merge) <- "integer"
    n1 <- nrow(merge)
    n <- n1 + 1L
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
    graphics:::plotHclust(n1, merge, height, order(x$order), 
        hang, labels, ...)
    if (axes) 
        axis(2, at = pretty(range(height)), ...)
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
