"summary.humpfit" <-
    function(object, ...)
{
    p <- coef(object)
    se <- sqrt(diag(solve(object$nlm$hessian)))
    est <- cbind(p, se)
    colnames(est) <- c("Estimate", "Std. Error")
    covmat <- solve(object$nlm$hessian)
    dg <- sqrt(diag(covmat))
    cormat <- covmat/outer(dg, dg)
    colnames(cormat) <- names(p)
    rownames(cormat) <- names(p)
    aic <- AIC(object)
    bic <- AIC(object, k = log(length(object$y)))
    out <- list(est = est, aic = aic, bic = bic, family = family(object)$family,
                deviance = deviance(object), df.residual = df.residual(object),
                correlation = cormat, iter = object$nlm$iterations,
                code = object$nlm$code)
    class(out) <- "summary.humpfit"
    out
}
