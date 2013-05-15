plot.JMbayes <-
function (x, which = c("trace", "autocorr"), param = c("betas", "tau", "inv.D", 
        "gammas", "alphas", "Dalphas", "sigma.t", "Bs.gammas"), ask = TRUE, ...) {
    if (!inherits(x, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    which <- match.arg(which)
    param <- match.arg(param, several.ok = TRUE)
    codaList <- x$codaFit
    codaList[] <- lapply(codaList, function (k) {
        ind <- sapply(paste("^", param, sep = ""), grep, x = colnames(k))
        if (!is.null(ind$Bs.gammas) && !is.null(x$coefficients$gammas))
            ind$gammas <- setdiff(ind$gammas, ind$Bs.gammas)
        k[, unlist(ind)]
    })
    if (which == "trace") {
        plot(codaList, ask = ask)
    } else {
        autocorr.plot(codaList, ask = ask)
    }
    invisible()
}
