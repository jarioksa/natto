`distconstrain` <-
    function(formula, data, residuals = FALSE, squared = FALSE)
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
    if (squared && any(dis < 0))
        warning("some dissimilarities were negative, use 'squared = TRUE'?")
    if (!squared)
        dis <- sqrt(dis)
    dis
}
