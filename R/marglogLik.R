marglogLik <-
function (object, newdata, idVar = "id", method = "BFGS", control = NULL) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    survMod <- object$survMod
    timeVar <- object$timeVar
    interFact <- object$interFact
    robust <- object$robust
    df <- object$df
    param <- object$param
    extraForm <- object$extraForm
    indFixed <- extraForm$indFixed
    indRandom <- extraForm$indRandom
    TermsX <- object$termsYx
    TermsZ <- object$termsYz
    TermsX.deriv <- object$termsYx.deriv
    TermsZ.deriv <- object$termsYz.deriv
    mfX <- model.frame(TermsX, data = newdata)
    mfZ <- model.frame(TermsZ, data = newdata)
    formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
    formYz <- object$formYz
    na.ind <- as.vector(attr(mfX, "na.action"))
    na.ind <- if (is.null(na.ind)) {
        rep(TRUE, nrow(newdata))
    } else {
        !seq_len(nrow(newdata)) %in% na.ind
    }
    id <- as.numeric(unclass(newdata[[idVar]]))
    id <- id. <- match(id, unique(id))
    id <- id[na.ind]
    y <- model.response(mfX)
    X <- model.matrix(formYx, mfX)
    Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
    TermsT <- object$termsT
    data.id <- newdata[!duplicated(id), ]
    idT <- data.id[[idVar]]
    idT <- match(idT, unique(idT))
    mfT <- model.frame(delete.response(TermsT), data = data.id)
    formT <- if (!is.null(kk <- attr(TermsT, "specials")$strata)) {
        strt <- eval(attr(TermsT, "variables"), data.id)[[kk]]
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
        tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
        reformulate(attr(tt, "term.labels"))
    } else {
        tt <- attr(delete.response(TermsT), "term.labels")
        if (length(tt)) reformulate(tt) else reformulate("1")
    }
    W <- model.matrix(formT, mfT)
    WintF.vl <- WintF.sl <- as.matrix(rep(1, nrow(data.id)))
    if (!is.null(interFact)) {
        if (!is.null(interFact$value))
            WintF.vl <- model.matrix(interFact$value, data = data.id)
        if (!is.null(interFact$slope))
            WintF.sl <- model.matrix(interFact$slope, data = data.id)
    }
    last.time <- tapply(newdata[[timeVar]], id., tail, n = 1)
    n.tp <- length(last.time)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    if (survMod == "spline-PH") {
        if (ncww == 1) {
            W <- NULL
            ncww <- 0
        } else {
            W <- W[, -1, drop = FALSE]
            ncww <- ncww - 1
        }
    }
    lag <- object$y$lag
    # list of parameters
    list.thetas <- c(object$modes, list(ranef = rbind(rep(0, ncz)) ))
    list.thetas$sigma <- log(list.thetas$sigma)
    list.thetas$D <- chol.transf(list.thetas$D)
    if (!is.null(list.thetas$sigma.t))
        list.thetas$sigma.t <- log(list.thetas$sigma.t)
    thetas.b <- unlist(as.relistable(list.thetas))
    # construct model matrices to calculate the log.posterior
    environment(ModelMats) <- environment()
    survMats.last <- vector("list", n.tp)
    for (i in seq_along(last.time)) {
        survMats.last[[i]] <- ModelMats(last.time[i], ii = i)
    }
    logPost <- function (thetas.b, ii, transform = TRUE) {
        # parameters
        tht <- relist(thetas.b, list.thetas)
        betas <- tht$betas        
        gammas <- tht$gammas
        alphas <- tht$alphas
        Dalphas <- tht$Dalphas
        Bs.gammas <- tht$Bs.gammas
        b <- tht$ranef
        if (transform) {
            sigma <- exp(tht$sigma)
            D <- chol.transf(tht$D)
            sigma.t <- if (is.null(tht$sigma.t)) NULL else exp(tht$sigma.t)
        } else {
            sigma <- tht$sigma
            D <- matrix(0, ncz, ncz)
            D[lower.tri(D, TRUE)] <- tht$D
            D <- D + t(D)
            diag(D) <- diag(D) / 2
            sigma.t <- tht$sigma.t
        }
        # log-likelihood contributions
        id.i <- id %in% ii
        idT.i <- idT %in% ii
        X.i <- X[id.i, , drop = FALSE]
        Z.i <- Z[id.i, , drop = FALSE]
        mu.y <- as.vector(X.i %*% betas) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
        logY <- if (!robust) dnorm(y[id.i], mu.y, sigma, TRUE) else dgt(y[id.i], mu.y, sigma, df, TRUE)
        log.p.yb <- sum(logY)
        log.p.b <- dmvnorm(b, rep(0, ncol(Z)), D, TRUE)
        st <- survMats.last[[ii]]$st
        wk <- survMats.last[[ii]]$wk
        P <- survMats.last[[ii]]$P
        Xs <- survMats.last[[ii]]$Xs
        Zs <- survMats.last[[ii]]$Zs
        Xs.deriv <- survMats.last[[ii]]$Xs.deriv
        Zs.deriv <- survMats.last[[ii]]$Zs.deriv
        Ws.intF.vl <- survMats.last[[ii]]$Ws.intF.vl
        Ws.intF.sl <- survMats.last[[ii]]$Ws.intF.sl
        ind <- survMats.last[[ii]]$ind
        if (param %in% c("td-value", "td-both"))
            Ys <- as.vector(Xs %*% betas + rowSums(Zs * rep(b, each = nrow(Zs))))
        if (param %in% c("td-extra", "td-both"))
            Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
                rowSums(Zs.deriv * rep(b[indRandom], each = nrow(Zs)))
        tt <- switch(param,
            "td-value" = c(Ws.intF.vl %*% alphas) * Ys, 
            "td-extra" = c(Ws.intF.sl %*% Dalphas) * Ys.deriv,
            "td-both" = c(Ws.intF.vl %*% alphas) * Ys + c(Ws.intF.sl %*% Dalphas) * Ys.deriv,
            "shared-RE" = rep(sum(b * alphas), length(st)))
        eta.tw <- if (!is.null(W)) {
                as.vector(W[ii, , drop = FALSE] %*% gammas)
        } else 0
        log.survival <- if (survMod == "weibull-PH") {
            Vi <- exp(log(sigma.t) + (sigma.t - 1) * log(st) + tt)
            - exp(eta.tw) * P * sum(wk * Vi)
        } else if (survMod == "spline-PH") {
            kn <- object$control$knots
            W2s <- splineDesign(unlist(kn, use.names = FALSE), st, 
                ord = object$control$ordSpline, outer.ok = TRUE)
            Vi <- exp(c(W2s %*% Bs.gammas) + tt)
            idT <- rep(seq_along(P), each = 15)
            - sum(exp(eta.tw) * P * tapply(wk * Vi, idT, sum))
        }
        if (all(st == 0))
            log.survival <- 1
        logLik <- log.p.yb + log.survival + log.p.b
        # Priors
        priors <- object$priors
        log.betas <- dmvnorm(betas, priors$priorMean.betas, 
            solve(priors$priorTau.betas), log = TRUE)
        log.tau <- dgamma(1 / (sigma^2), priors$priorA.tau, 
            priors$priorB.tau, log = TRUE)
        log.D <- dwish(D, priors$priorR.D, priors$priorK.D, log = TRUE)
        logPrior <- log.betas + log.tau + log.D
        if (!is.null(gammas)) {
            ind <- colSums(object$x$W == 0) == nrow(object$x$W)
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
        - as.vector(logLik + logPrior)
    }
    score.logPost <- function (thetas.b, ii, transform = TRUE) {
        fd(thetas.b, logPost, ii = ii, transform = transform)
    }
    w <- numeric(length(last.time))
    con <- list(maxit = 200, parscale = rep(0.01, length(thetas.b)))
    con[names(control)] <- control
    for (i in seq_along(last.time)) {
        w[i] <- tryCatch({
            opt <- optim(thetas.b, logPost, score.logPost, ii = i, transform = TRUE, 
                method = method, control = con, hessian = TRUE)
            opt.thetas <- relist(opt$par, list.thetas)
            opt.thetas$sigma <- exp(opt.thetas$sigma)
            opt.thetas$D <- chol.transf(opt.thetas$D)
            opt.thetas$D <- opt.thetas$D[lower.tri(opt.thetas$D, TRUE)]
            opt.thetas$sigma.t <- if (is.null(opt.thetas$sigma.t)) NULL else exp(opt.thetas$sigma.t)
            opt.thetas <- unlist(as.relistable(opt.thetas))
            H <- opt$hessian #fd.vec(opt.thetas, score.logPost, ii = i, transform = FALSE, eps = 1e-06)
            as.vector(0.5 * length(unlist(list.thetas)) * log(2 * pi) - 
                0.5 * determinant(H)$modulus - opt$value)
        }, error = function (e) NA)
    }
    w
}
