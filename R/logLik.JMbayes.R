logLik.JMbayes <-
function (object, thetas, b, priors = TRUE, marginal.b = TRUE, marginal.thetas = FALSE, 
        full.Laplace = FALSE, useModes = TRUE, ...) {
    if (!inherits(object, "JMbayes"))
        stop("'object' must inherit from class JMbayes.")
    if (missing(thetas))
        thetas <- object$coefficients
    if (missing(b))
        b <- ranef(object)
    if (!is.list(thetas) || length(thetas) != length(object$coefficients))
        stop("'thetas' must be a list with the model's parameters with the same structure as ",
            "'object$coefficients'.")    
    if (!is.matrix(b) || dim(b) != dim(ranef(object)))
        stop("'b' must be a numeric matrix with random effects values with the same ",
            "dimensions as 'ranef(object)'.")
    # data and settings
    survMod <- object$survMod
    timeVar <- object$timeVar
    robust <- object$robust
    robust.b <- object$robust.b
    lag <- object$y$lag
    df <- object$df
    df.b <- object$df.b
    param <- object$param
    indFixed <- object$extraForm$indFixed
    indRandom <- object$extraForm$indRandom
    id <- object$id
    id <- match(id, unique(id))
    y <- object$y$y
    logT <- object$y$logT
    Time <- exp(logT)
    event <- object$y$event
    n <- length(logT)
    X <- object$x$X
    Z <- object$x$Z
    W <- object$x$W
    if (any(ind <- colSums(W == 0) == nrow(W)))
        W <- W[, !ind, drop = FALSE]
    W2 <- object$x$W2
    W2s <- object$x$W2s
    Xtime <- object$x$Xtime
    Ztime <- object$x$Ztime
    Xtime.deriv <- object$x$Xtime.deriv
    Ztime.deriv <- object$x$Ztime.deriv
    Xs <- object$x$Xs
    Zs <- object$x$Zs
    Xs.deriv <- object$x$Xs.deriv
    Zs.deriv <- object$x$Zs.deriv
    st <- object$x$st
    wk <- object$x$wk
    P <- object$x$P
    # parameters
    betas <- thetas$betas
    sigma <- thetas$sigma
    D <- thetas$D
    gammas <- thetas$gammas
    alphas <- thetas$alphas
    Dalphas <- thetas$Dalphas
    sigma.t <- thetas$sigma.t
    Bs.gammas <- thetas$Bs.gammas
    # log-likelihood
    h <- function (b, individuals = NULL) {
        mu.y <- c(X %*% betas) + rowSums(Z * b[id, , drop = FALSE])
        log.p.y.b <- if (!robust) dnorm(y, mu.y, sigma, log = TRUE) else dgt(y, mu.y, sigma, df, log = TRUE)
        log.p.y.b <- tapply(log.p.y.b, id, sum)
        eta.t <- if (!is.null(gammas)) c(W %*% gammas) else rep(0, n)
        id.GK <- rep(seq_along(logT), each = 15)
        wk.long <- rep(wk, n)
        if (param %in% c("td-value", "td-both")) {
            Y <- c(Xtime %*% betas) + rowSums(Ztime * b)
            Ys <- c(Xs %*% betas) + rowSums(Zs * b[id.GK, , drop = FALSE])
        }
        if (param %in% c("td-extra", "td-both")) {
            Yderiv <- c(Xtime.deriv %*% betas[indFixed]) + rowSums(Ztime.deriv * b[, indRandom, drop = FALSE])
            Ys.deriv <- c(Xs.deriv %*% betas[indFixed]) + rowSums(Zs.deriv * b[id.GK, indRandom, drop = FALSE])
        }
        longSurv <- switch(param,
            "td-value" = alphas * Y, 
            "td-extra" = Dalphas * Yderiv,
            "td-both" = alphas * Y + Dalphas * Yderiv,
            "shared-RE" = c(b %*% alphas))
        longSurv.s <- switch(param,
            "td-value" = alphas * Ys, 
            "td-extra" = Dalphas * Ys.deriv,
            "td-both" = alphas * Ys + Dalphas * Ys.deriv,
            "shared-RE" = c(b %*% alphas)[id.GK])
        if (survMod == "weibull-PH") {
            log.hazard <- log(sigma.t) + (sigma.t - 1) * logT + eta.t + longSurv
            log.survival <- - exp(eta.t) * P * tapply(wk.long * exp(log(sigma.t) + 
                (sigma.t - 1) * log(c(t(st))) + longSurv.s), id.GK, sum)
        } else {
            log.hazard <- c(W2 %*% Bs.gammas) + eta.t + longSurv
            log.survival <- - exp(eta.t) * P * tapply(wk.long * exp(c(W2s %*% Bs.gammas) + 
                longSurv.s), id.GK, sum)
        }
        log.p.t.b <- event * log.hazard + log.survival
        if (!is.null(individuals))
            log.p.y.b[individuals] + log.p.t.b[individuals]
        else
            log.p.y.b + log.p.t.b
    }
    logLik <- if (!marginal.b) {
        log.p.b <- if (!robust.b) dmvnorm(b, rep(0, ncol(b)), D, log = TRUE)
            else dmvt(b, rep(0, ncol(b)), D, df.b, log = TRUE)
        sum(h(b) + log.p.b, na.rm = TRUE)
    } else {
        mean.b <- ranef(object)
        var.b <- attr(ranef(object, postVar = TRUE), "postVar")
        if (full.Laplace) {
            optFun <- function (b, id) {
                b. <- ranef(object)
                b.[id, ] <- b
                log.p.b <- if (!robust.b) dmvnorm(b., rep(0, ncol(b.)), D, log = TRUE)
                    else dmvt(b., rep(0, ncol(b.)), D, df.b, log = TRUE)
                - h(b., id) - log.p.b[id] 
            }
            for (i in seq_len(n)) {
                opt <- optim(mean.b[i, ], optFun, id = i, method = "BFGS", 
                    hessian = TRUE)
                mean.b[i, ] <- opt$par
                var.b[[i]] <- solve(opt$hessian)
            }
        }
        log.p.b <- if (!robust.b) dmvnorm(mean.b, rep(0, ncol(b)), D, log = TRUE)
            else dmvt(mean.b, rep(0, ncol(b)), D, df.b, log = TRUE)
        log.dets.var.b <- sapply(var.b, function (x) determinant(x)$modulus)
        sum(0.5 * ncol(b) * log(2 * pi) + 0.5 * log.dets.var.b + 
            h(mean.b) + log.p.b, na.rm = TRUE)
    }
    # priors
    if (priors) {
        priors <- object$priors
        log.betas <- dmvnorm(betas, priors$priorMean.betas, 
            solve(priors$priorTau.betas), log = TRUE)
        log.tau <- dgamma(1 / (sigma^2), priors$priorA.tau, 
            priors$priorB.tau, log = TRUE)
        log.D <- dwish(D, priors$priorR.D, priors$priorK.D, log = TRUE)
        logPrior <- log.betas + log.tau + log.D
        if (!is.null(gammas)) {
            log.gammas <- dmvnorm(gammas, priors$priorMean.gammas[!ind], 
                solve(priors$priorTau.gammas)[!ind, !ind], log = TRUE)
            logPrior <- logPrior + log.gammas
        }
        if (!is.null(alphas)) {
            log.alphas <- dmvnorm(alphas, priors$priorMean.alphas, 
                solve(priors$priorTau.alphas), log = TRUE)
            logPrior <- logPrior + log.alphas
        }
        if (!is.null(Dalphas)) {
            log.Dalphas <- dmvnorm(Dalphas, priors$priorMean.Dalphas, 
                solve(priors$priorTau.Dalphas), log = TRUE)
            logPrior <- logPrior + log.Dalphas
        }
        if (!is.null(sigma.t)) {
            log.nu <- dgamma(sigma.t, priors$priorA.sigma.t, priors$priorB.sigma.t, log = TRUE)
            logPrior <- logPrior + log.nu
        }
        if (!is.null(Bs.gammas)) {
            log.Bs.gammas <- dmvnorm(Bs.gammas, priors$priorMean.Bs.gammas, 
                solve(priors$priorTau.Bs.gammas), log = TRUE)
            logPrior <- logPrior + log.Bs.gammas
        }
        logLik <- logLik + logPrior
    }
    if (marginal.thetas) {
        tht <- if (useModes) object$modes else object$coefficients
        lL <- logLik(object, thetas = tht)
        var.thetas <- hessian.JMbayes(object, thetas = tht)
        tht$D <- tht$D[lower.tri(tht$D, TRUE)]
        nthetas <- length(unlist(tht))
        logLik <- 0.5 * nthetas * log(2 * pi) - 
            0.5 * determinant(var.thetas)$modulus + lL
    }
    out <- as.vector(logLik)
    attr(out, "df") <- nrow(object$vcov)
    attr(out, "nobs") <- n
    class(out) <- "logLik"
    out
}
