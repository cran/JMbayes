plot.JMbayes <-
function (x, which = c("trace", "autocorr", "CPO"), 
                          param = c("betas", "sigma", "D", "gammas", "alphas", "Dalphas", 
                                    "shapes", "Bs.gammas", "tauBs"), 
                          ask = TRUE, ...) {
    if (!inherits(x, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    which <- match.arg(which)
    if (which %in% c("trace", "autocorr")) {
        param <- match.arg(param, several.ok = TRUE)
        if (any(param == "D")) {
            keepD <- c(lower.tri(x$postMeans$D, TRUE))
            x$mcmc$D <- x$mcmc$D[, keepD]
        }
        pp <- do.call(cbind, x$mcmc[param])
        nams <- colnames(pp)
        op <- if (ask) par(mfrow = c(2, 2), ask = ask) else par(mfrow = c(4, 2))
        if (which == "trace") {   
            for (i in 1:ncol(pp))
                plot(pp[, i], type = "l", xlab = "iterations", ylab = nams[i])
        } else {
            for (i in 1:ncol(pp))
                acf(pp[, i], ylab = nams[i], main = paste("Series", nams[i]))
        }
        par(op)        
    } else {
        n <- length(x$CPO)
        matplot(matrix(1:n, 2, n, TRUE), rbind(rep(0, n), x$CPO), type = "l",
                xlab = "Subjects", ylab = "CPO", main = deparse(substitute(x)), ...)
    }
    invisible()
}
