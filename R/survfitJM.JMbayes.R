survfitJM.JMbayes <-
function (object, newdata, idVar = "id", simulate = TRUE, survTimes = NULL, 
            last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6, weight = rep(1, nrow(newdata)), ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    if (is.null(survTimes) || !is.numeric(survTimes))
        survTimes <- seq(min(exp(object$y$logT)), 
            max(exp(object$y$logT)) + 0.1, length.out = 35)
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
    obs.times <- split(newdata[[timeVar]][na.ind], id)
    last.time <- if (is.null(last.time)) {
        tapply(newdata[[timeVar]], id., tail, n = 1)
    } else if (is.character(last.time) && length(last.time) == 1) {
        tapply(newdata[[last.time]], id., tail, n = 1)
    } else if (is.numeric(last.time) && length(last.time) == nrow(data.id)) {
        last.time
    } else {
        stop("\nnot appropriate value for 'last.time' argument.")
    }
    times.to.pred <- lapply(last.time, function (t) survTimes[survTimes > t])
    n <- object$n
    n.tp <- length(last.time)
    ncx <- ncol(X)
    ncz <- ncol(Z)
    ncww <- ncol(W)
    lag <- object$y$lag
    betas <- object$coefficients$betas
    sigma <- object$coefficients$sigma
    D <- object$coefficients$D
    gammas <- object$coefficients$gammas
    alpha <- object$coefficients$alpha
    Dalpha <- object$coefficients$Dalpha
    sigma.t <- object$coefficients$sigma.t
    Bs.gammas <- object$coefficients$Bs.gammas
    list.thetas <- list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha, 
            Dalpha = Dalpha, sigma.t = sigma.t, Bs.gammas = Bs.gammas, D = D)
    if (survMod == "spline-PH") {
        if (ncww == 1) {
            W <- NULL
            ncww <- 0
        } else {
            W <- W[, -1, drop = FALSE]
            ncww <- ncww - 1
        }
    }
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    Var.thetas <- vcov(object)
    environment(log.posterior.b) <- environment(S.b) <- environment(ModelMats) <- environment()
    # construct model matrices to calculate the survival functions
    obs.times.surv <- split(data.id[[timeVar]], idT)
    survMats <- survMats.last <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        survMats[[i]] <- lapply(times.to.pred[[i]], ModelMats, ii = i)
        survMats.last[[i]] <- ModelMats(last.time[i], ii = i)
    }
    # calculate the Empirical Bayes estimates and their (scaled) variance
    modes.b <- matrix(0, n.tp, ncz)
    Vars.b <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
        betas.new <- betas
        sigma.new <- sigma
        D.new <- D
        gammas.new <- gammas
        alpha.new <- alpha
        Dalpha.new <- Dalpha
        sigma.t.new <- sigma.t
        Bs.gammas.new <- Bs.gammas
        ff <- function (b, y, tt, mm, i) -log.posterior.b(b, y, Mats = tt, survMod = mm, ii = i)
        opt <- try(optim(rep(0, ncz), ff, y = y, tt = survMats.last, mm = survMod, i = i, 
            method = "BFGS", hessian = TRUE), TRUE)
        if (inherits(opt, "try-error")) {
            gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, mm = mm, i = i)
            opt <- optim(rep(0, ncz), ff, gg, y = y, tt = survMats.last, mm = survMod, 
                i = i, method = "BFGS", hessian = TRUE, control = list(parscale = rep(0.1, ncz)))
        } 
        modes.b[i, ] <- opt$par
        Vars.b[[i]] <- scale * solve(opt$hessian)
    }
    if (!simulate) {
        res <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            S.last <- S.b(last.time[i], modes.b[i, ], i, survMats.last[[i]])
            S.pred <- numeric(length(times.to.pred[[i]]))
            for (l in seq_along(S.pred))
                S.pred[l] <- S.b(times.to.pred[[i]][l], modes.b[i, ], i, survMats[[i]][[l]])
            res[[i]] <- cbind(times = times.to.pred[[i]], predSurv = weight[i] * S.pred / S.last)
            rownames(res[[i]]) <- seq_along(S.pred) 
        }
    } else {
        out <- vector("list", M)
        success.rate <- matrix(FALSE, M, n.tp)
        b.old <- b.new <- modes.b
        if (n.tp == 1)
            dim(b.old) <- dim(b.new) <- c(1, ncz)
        codaFit <- do.call("rbind", object$codaFit)
        if (M > nrow(codaFit)) {
            warning("'M' cannot be set greater than ", nrow(codaFit))
            M <- nrow(codaFit)
            out <- vector("list", M)
            success.rate <- matrix(FALSE, M, n.tp)
        }
        codaFit <- codaFit[sample(nrow(codaFit), M), , drop = FALSE]
        codaFit <- cbind(codaFit, "sigma" = 1 / codaFit[, "tau"])
        ind.thetas <- sapply(names(list.thetas), grep, x = colnames(codaFit), fixed = TRUE)
        if (!is.null(ind.thetas$Bs.gammas) && !is.null(object$coefficients$gammas))
            ind.thetas$gammas <- setdiff(ind.thetas$gammas, ind.thetas$Bs.gammas)
        if (param %in% c("td-extra", "td-both"))
            ind.thetas$D <- setdiff(ind.thetas$D, ind.thetas$Dalpha)
        if (param == "td-both")
            ind.thetas$alpha <- setdiff(ind.thetas$alpha, ind.thetas$Dalpha)
        if (param == "td-extra")
            ind.thetas$alpha <- NULL
        if (param %in% c("td-value", "shared-RE"))
            ind.thetas$Dalpha <- NULL
        if (survMod == "weibull-PH") {
            ind.thetas$sigma <- setdiff(ind.thetas$sigma, ind.thetas$sigma.t)
            ind.thetas$Bs.gammas <- NULL
        }
        for (m in 1:M) {
            # Step 1: extract parameter values
            thetas.new <- lapply(ind.thetas, function (k) codaFit[m, k])
            thetas.new$D <- solve(matrix(thetas.new$D, ncz, ncz, TRUE))
            betas.new <- thetas.new$betas
            sigma.new <- thetas.new$sigma
            gammas.new <- thetas.new$gammas
            alpha.new <- thetas.new$alpha
            Dalpha.new <- thetas.new$Dalpha
            D.new <- thetas.new$D
            sigma.t.new <- thetas.new$sigma.t
            Bs.gammas.new <- thetas.new$Bs.gammas
            SS <- vector("list", n.tp)
            for (i in seq_len(n.tp)) {
                # Step 2: simulate new random effects values
                proposed.b <- rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
                dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
                dmvt.proposed <- dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
                a <- min(exp(log.posterior.b(proposed.b, y, survMats.last, survMod, ii = i) + dmvt.old - 
                        log.posterior.b(b.old[i, ], y, survMats.last, survMod, ii = i) - dmvt.proposed), 1)
                ind <- runif(1) <= a
                success.rate[m, i] <- ind
                if (!is.na(ind) && ind)
                    b.new[i, ] <- proposed.b
                # Step 3: compute Pr(T > t_k | T > t_{k - 1}; theta.new, b.new)
                S.last <- S.b(last.time[i], b.new[i, ], i, survMats.last[[i]])
                S.pred <- numeric(length(times.to.pred[[i]]))
                for (l in seq_along(S.pred))
                    S.pred[l] <- S.b(times.to.pred[[i]][l], b.new[i, ], i, survMats[[i]][[l]])
                SS[[i]] <- weight[i] * S.pred / S.last
            }
            b.old <- b.new
            out[[m]] <- SS
        }
        res <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            rr <- sapply(out, "[[", i)
            if (!is.matrix(rr))
                rr <- rbind(rr)
            res[[i]] <- cbind(
                times = times.to.pred[[i]],
                "Mean" = rowMeans(rr, na.rm = TRUE),
                "Median" = apply(rr, 1, median, na.rm = TRUE),
                "Lower" = apply(rr, 1, quantile, probs = CI.levels[1], na.rm = TRUE),
                "Upper" = apply(rr, 1, quantile, probs = CI.levels[2], na.rm = TRUE)
            )
            rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
        }
    }
    y <- split(y, id)
    newdata. <- do.call(rbind, mapply(function (d, t) {
        d. <- rbind(d, d[nrow(d), ])
        d.[[timeVar]][nrow(d.)] <- t
        d.
    }, split(newdata, id.), last.time, SIMPLIFY = FALSE))
    id. <- as.numeric(unclass(newdata.[[idVar]]))
    id. <- match(id., unique(id.))
    mfX. <- model.frame(delete.response(TermsX), data = newdata.)
    mfZ. <- model.frame(TermsZ, data = newdata.)
    X. <- model.matrix(formYx, mfX.)
    Z. <- model.matrix(formYz, mfZ.)
    fitted.y <- split(c(X. %*% betas) + rowSums(Z. * modes.b[id., , drop = FALSE]), id.)
    names(res) <- names(y) <- names(last.time) <- names(obs.times) <- unique(unclass(newdata[[idVar]]))
    res <- list(summaries = res, survTimes = survTimes, last.time = last.time, 
        obs.times = obs.times, y = y, 
        fitted.times = split(newdata.[[timeVar]], factor(newdata.[[idVar]])), 
        fitted.y = fitted.y, ry = range(object$y$y, na.rm = TRUE))
    if (simulate) {
        res$full.results <- out
        res$success.rate <- success.rate
    }
    class(res) <- "survfit.JMbayes"
    res
}
