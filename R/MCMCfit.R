MCMCfit <-
function (y, x, param, extraForm, baseHaz, estimateWeightFun, initials, priors, 
                     scales, Funs, Covs, Data, control, df.RE) {
    # extract data longitudinal
    y.long <- dropAttr(y$y)
    X <- dropAttr(x$X)
    Z <- dropAttr(x$Z)
    # extract data survival
    Time <- dropAttr(y$Time)
    event <- dropAttr(y$event)
    W <- dropAttr(x$W)
    notNullW <- !is.null(W)
    Ws <- dropAttr(x$Ws)
    W2 <- dropAttr(x$W2)
    W2s <- dropAttr(x$W2s)
    Xtime <- dropAttr(x$Xtime)
    Ztime <- dropAttr(x$Ztime)
    Xtime.extra <- dropAttr(x$Xtime.extra)
    Ztime.extra <- dropAttr(x$Ztime.extra)
    Xs <- dropAttr(x$Xs)
    Zs <- dropAttr(x$Zs)
    Xs.extra <- dropAttr(x$Xs.extra)
    Zs.extra <- dropAttr(x$Zs.extra)
    Xu <- dropAttr(x$Xu)
    Zu <- dropAttr(x$Zu)
    # extract indices
    n <- length(Time)
    ns <- nrow(W2s)
    nu <- nrow(Xu)
    ncZ <- ncol(Z)
    nrZ <- nrow(Z)
    nrZtime <- nrow(Ztime)
    nrZs <- nrow(Zs)
    ncZ.extra <- ncol(Ztime.extra)
    nrZtime.extra <- nrow(Ztime.extra)
    nrZs.extra <- nrow(Zs.extra)
    nrZu <- nrow(Zu)
    id <- dropAttr(y$id)
    id.GK <- dropAttr(y$id.GK)
    indBetas <- dropAttr(y$indBetas)
    indBetas2 <- rep(indBetas, each = n)
    iF <- dropAttr(extraForm$indFixed)
    iR <- dropAttr(extraForm$indRandom)
    lag <- y$lag
    w <- rep(dropAttr(x$wk), n)
    P <- dropAttr(x$P)
    st <- c(t(dropAttr(x$st)))
    if (estimateWeightFun) {
        nshapes <- length(initials$shapes)
        seq.nshapes <- seq_len(nshapes)
        weightFun <- Funs$weightFun
        id.GK2 <- dropAttr(y$id.GK2)
        id.GKu <- rep(id.GK, each = length(x$wk))
        w2 <- rep(dropAttr(x$wk), ns)
        P2 <- dropAttr(x$P2)
        st2 <- c(t(dropAttr(x$st2)))
        max.time <- max(Time)
        u.idGK <- Time[id.GK] - st
        u.idGK2 <- st[id.GK2] - st2
    }
    paramValue <- (param %in% c("td-value", "td-both")) && !estimateWeightFun
    paramExtra <- param %in% c("td-extra", "td-both")
    paramRE <- param %in% c("shared-betasRE", "shared-RE")
    paramSharedRE <- param == "shared-RE"
    baseHazP <- baseHaz == "P-splines"
    paramValueRE <- (paramValue || paramRE)
    estimateAlphas <- paramValueRE || estimateWeightFun
    notestimateWeightFun <- !estimateWeightFun
    # extract initial values
    init.betas <- betas <- dropAttr(initials$betas)
    init.tau <- tau <- dropAttr(initials$tau)
    init.b <- b <- dropAttr(initials$b)
    init.invD <- invD <- dropAttr(initials$invD)
    init.gammas <- gammas <- dropAttr(initials$gammas)
    init.Bs.gammas <- Bs.gammas <- dropAttr(initials$Bs.gammas)
    init.tauBs <- tauBs <- dropAttr(initials$tauBs)
    init.alphas <- alphas <- dropAttr(initials$alphas)
    init.Dalphas <- Dalphas <- dropAttr(initials$Dalphas)
    init.shapes <- shapes <- dropAttr(initials$shapes)
    # dimensions of parameters
    nbetas <- length(betas)
    nRE <- rep(ncZ, n)
    ngammas <- length(gammas)
    nBs.gammas <- length(Bs.gammas)
    nalphas <- length(alphas)
    nDalphas <- length(Dalphas)
    # extract Funs
    densLong <- Funs$densLong
    hasScale <- Funs$hasScale
    densRE <- Funs$densRE
    transFun.value <- Funs$transFun.value
    transFun.extra <- Funs$transFun.extra
    # Data sets
    data <- Data$data
    data.id <- Data$data.id
    data.s <- Data$data.s
    data.u <- Data$data.u
    # define priors
    priorMean.betas <- priors$priorMean.betas
    priorTau.betas <- priors$priorTau.betas
    log.prior.betas <- function (betas) {
        dmvnorm(betas, priorMean.betas, invSigma = priorTau.betas, log = TRUE)
    }
    priorA.tau <- priors$priorA.tau
    priorB.tau <- priors$priorB.tau
    log.prior.tau <- function (tau) {
        dgamma(tau, priorA.tau, priorB.tau)
    }
    priorR.invD <- priors$priorR.invD
    priorK.invD <- priors$priorK.invD
    log.prior.invD <- function (invD) {
        dwish(invD, priorR.invD, priorK.invD, log = TRUE)
    }
    priorMean.gammas <- priors$priorMean.gammas
    priorTau.gammas <- priors$priorTau.gammas
    log.prior.gammas <- function (gammas) {
        dmvnorm(gammas, priors$priorMean.gammas, invSigma = priorTau.gammas, log = TRUE)
    }
    priorMean.Bs.gammas <- priors$priorMean.Bs.gammas
    priorTau.Bs.gammas <- priors$priorTau.Bs.gammas
    log.prior.Bs.gammas <- function (Bs.gammas) {
        if (baseHazP)
            priorTau.Bs.gammas <- tauBs * priorTau.Bs.gammas
        dmvnorm(Bs.gammas, priorMean.Bs.gammas, invSigma = priorTau.Bs.gammas, log = TRUE)
    }
    priorA.tauBs <- priors$priorA.tauBs
    priorB.tauBs <- priors$priorB.tauBs
    priorMean.alphas <- priors$priorMean.alphas
    priorTau.alphas <- priors$priorTau.alphas
    log.prior.alphas <- function (alphas) {
        dmvnorm(alphas, priorMean.alphas, invSigma = priorTau.alphas, log = TRUE)
    }
    priorMean.Dalphas <- priors$priorMean.Dalphas
    priorTau.Dalphas <- priors$priorTau.Dalphas
    log.prior.Dalphas <- function (Dalphas) {
        dmvnorm(Dalphas, priorMean.Dalphas, invSigma = priorTau.Dalphas, log = TRUE)
    }
    priorshape1Fun <- control$priorShapes$shape1
    priorshape1.low <- priors$priorshape1[1L]
    priorshape1.upp <- priors$priorshape1[2L]
    log.prior.shape1 <- function (shape1) {
        priorshape1Fun(shape1, priorshape1.low, priorshape1.upp, log = TRUE)
    }
    priorshape2Fun <- control$priorShapes$shape2
    priorshape2.low <- priors$priorshape2[1L]
    priorshape2.upp <- priors$priorshape2[2L]
    log.prior.shape2 <- function (shape2) {
        priorshape2Fun(shape2, priorshape2.low, priorshape2.upp, log = TRUE)
    }
    priorshape3Fun <- control$priorShapes$shape3
    priorshape3.low <- priors$priorshape3[1L]
    priorshape3.upp <- priors$priorshape3[2L]
    log.prior.shape3 <- function (shape3) {
        priorshape3Fun(shape3, priorshape3.low, priorshape3.upp, log = TRUE)
    }
    # define posteriors
    logPost.betas <- function (betas){
        Xbetas <- drop(X %*% betas)
        eta.y <- Xbetas + Zb
        log.pyb <- fastSumID(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), id)
        log.prior <- log.prior.betas(betas)
        if (!paramRE) {
            Mtime <- numeric(n)
            Ms <- numeric(ns)
            if (paramValue) {
                Xtimebetas <- drop(Xtime %*% betas)
                Xsbetas <- drop(Xs %*% betas)
                vl <- transFun.value(Xtimebetas + Ztimeb, data.id)
                vls <- transFun.value(Xsbetas + Zsb, data.s)
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            if (paramExtra) {
                Xtime.extrabetas <- drop(Xtime.extra %*% betas[iF])
                Xs.extrabetas <- drop(Xs.extra %*% betas[iF])
                ex <- transFun.extra(Xtime.extrabetas + Ztime.extrab, data.id)
                exs <- transFun.extra(Xs.extrabetas + Zs.extrab, data.s)
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
                Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
            }
            if (estimateWeightFun) {
                Xsbetas <- drop(Xs %*% betas)
                vl <- transFun.value(P * fastSumID(wFun * (Xsbetas + Zsb), id.GK), data.id)
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Xubetas <- drop(Xu %*% betas)
                vls <- transFun.value(P2 * fastSumID(wFun2 * (Xubetas + Zub), id.GK2), data.s)
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            log.Surv <- Int <- P * fastSumID(w * exp(log.h0s + Ms), id.GK)
            if (notNullW)
                log.Surv <- expWgammas * log.Surv
            log.ptb <- event * Mtime - log.Surv
            list(log.post = sum(log.pyb, log.ptb, na.rm = TRUE) + log.prior,
                 Xbetas = Xbetas,
                 Xtimebetas = if (paramValue) Xtimebetas,
                 Xsbetas = if (estimateAlphas) Xsbetas,
                 Xtime.extrabetas = if (paramExtra) Xtime.extrabetas,
                 Xs.extrabetas = if (paramExtra) Xs.extrabetas,
                 Xubetas = if (estimateWeightFun) Xubetas,
                 vl = if (estimateAlphas) vl, vls = if (estimateAlphas) vls,
                 ex = if (paramExtra) ex, exs = if (paramExtra) exs, 
                 Ms = Ms, Mtime = Mtime, log.Surv = log.Surv, Int = Int)
        } else {
            if (paramSharedRE) {
                list(log.post = sum(log.pyb, na.rm = TRUE) + log.prior, Xbetas = Xbetas)
            } else {
                Mtime <- drop((betas[indBetas2] + b) %*% alphas)
                log.Surv <- exp(Mtime) * Int
                if (notNullW)
                    log.Surv <- expWgammas * log.Surv
                log.ptb <- event * Mtime - log.Surv
                list(log.post = sum(log.pyb, log.ptb, na.rm = TRUE) + log.prior, Xbetas = Xbetas,
                     log.Surv = log.Surv)
            }
        }
    }
    logPost.betas2 <- function () {
        log.pyb <- fastSumID(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), id)
        log.prior <- log.prior.betas(betas)
        if (!paramRE) {
            Mtime <- numeric(n)
            if (paramValue) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            if (paramExtra) {
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            }
            if (estimateWeightFun) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            log.ptb <- event * Mtime - log.Surv
            sum(log.pyb, log.ptb, na.rm = TRUE) + log.prior
        } else {
            if (paramSharedRE) {
                sum(log.pyb, na.rm = TRUE) + log.prior
            } else {
                Mtime <- drop((betas[indBetas2] + b) %*% alphas)
                log.ptb <- event * Mtime - log.Surv
                sum(log.pyb, log.ptb, na.rm = TRUE) + log.prior
            }
        }
    }
    logPost.tau <- function (tau) {
        log.pyb <- densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data)
        log.prior <- log.prior.tau(tau)
        sum(log.pyb, na.rm = TRUE) + log.prior
    }
    logPost.RE <- function (b) {
        Zb <- .rowSums(Z * b[id, , drop = FALSE], nrZ, ncZ)
        eta.y <- Xbetas + Zb
        log.pyb <- fastSumID(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), id)
        log.prior <- densRE(b, invD = invD, log = TRUE)
        if (!paramRE) {
            Mtime <- numeric(n)
            Ms <- numeric(ns)
            if (paramValue) {
                Zsb <- .rowSums(Zs * b[id.GK, , drop = FALSE], nrZs, ncZ)
                Ztimeb <- .rowSums(Ztime * b, nrZtime, ncZ)
                vl <- transFun.value(Xtimebetas + Ztimeb, data.id)
                vls <- transFun.value(Xsbetas + Zsb, data.s)
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            if (paramExtra) {
                Ztime.extrab <- .rowSums(Ztime.extra * b[, iR, drop = FALSE], nrZtime.extra, ncZ.extra)
                Zs.extrab <- .rowSums(Zs.extra * b[id.GK, iR, drop = FALSE], nrZs.extra, ncZ.extra)
                ex <- transFun.extra(Xtime.extrabetas + Ztime.extrab, data.id)
                exs <- transFun.extra(Xs.extrabetas + Zs.extrab, data.s)
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
                Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
            }
            if (estimateWeightFun) {
                Zsb <- .rowSums(Zs * b[id.GK, , drop = FALSE], nrZs, ncZ)
                vl <- transFun.value(P * fastSumID(wFun * (Xsbetas + Zsb), id.GK), data.id)
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Zub <- .rowSums(Zu * b[id.GKu, , drop = FALSE], nrZu, ncZ)
                vls <- transFun.value(P2 * fastSumID(wFun2 * (Xubetas + Zub), id.GK2), data.s)
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            log.Surv <- Int <- P * fastSumID(w * exp(log.h0s + Ms), id.GK)
            if (notNullW)
                log.Surv <- expWgammas * log.Surv
            log.ptb <- event * Mtime - log.Surv
            list(log.post = log.pyb + log.ptb + log.prior,
                 Zb = Zb, Int = Int, Ms = Ms, eta.y = eta.y, log.Surv = log.Surv,
                 Ztimeb = if (paramValue) Ztimeb,
                 Zsb = if (estimateAlphas) Zsb,
                 Ztime.extrab = if (paramExtra) Ztime.extrab,
                 Zs.extrabetas = if (paramExtra) Zs.extrab,
                 Zub = if (estimateWeightFun) Zub,
                 vl = if (estimateAlphas) vl,
                 vls = if (estimateAlphas) vls,
                 ex = if (paramExtra) ex,
                 exs = if (paramExtra) exs)
        } else {
            Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
            log.Surv <- exp(Mtime) * Int
            if (notNullW)
                log.Surv <- expWgammas * log.Surv
            log.ptb <- event * Mtime - log.Surv
            list(log.post = log.pyb + log.ptb + log.prior,
                 Zb = Zb, eta.y = eta.y, log.Surv = log.Surv)
        }
    }
    logPost.RE2 <- function () {
        log.pyb <- fastSumID(densLong(y.long, Xbetas + Zb, 1/sqrt(tau), log = TRUE, data), id)
        log.prior <- densRE(b, invD = invD, log = TRUE)
        if (!paramRE) {
            Mtime <- numeric(n)
            if (paramValue) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            if (paramExtra) {
                Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            }
            if (estimateWeightFun) {
                Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            }
            log.ptb <- event * Mtime - log.Surv
        } else {
            Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
            log.ptb <- event * Mtime - log.Surv
        }
        log.pyb + log.ptb + log.prior
    }
    logPost.invD <- function (invD) {
        log.pb <- densRE(b, invD = invD, log = TRUE, prop = FALSE)
        log.prior <- log.prior.invD(invD)
        sum(log.pb, na.rm = TRUE) + log.prior
    }
    logPost.gammas <- function (gammas){
        Wgammas <- drop(W %*% gammas)
        log.Surv <- exp(Wgammas) * Int
        if (paramRE) {
            Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
            log.Surv <- exp(Mtime) * log.Surv
        }
        log.ptb <- event * Wgammas - log.Surv
        log.prior <- log.prior.gammas(gammas)
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, expWgammas = exp(Wgammas),
             log.Surv = log.Surv)
    }
    logPost.Bs.gammas <- function (Bs.gammas) {
        W2sBs.gammas <- drop(W2s %*% Bs.gammas)
        log.Surv <- Int <- P * if (!paramRE) {
            fastSumID(w * exp(W2sBs.gammas + Ms), id.GK)
        } else {
            fastSumID(w * exp(W2sBs.gammas), id.GK)
        }
        if (paramRE) {
            Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
            log.Surv <- exp(Mtime) * log.Surv
        }
        if (notNullW)
            log.Surv <- expWgammas * log.Surv
        log.ptb <- event * drop(W2 %*% Bs.gammas) - log.Surv
        log.prior <- log.prior.Bs.gammas(Bs.gammas)
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.h0s = W2sBs.gammas,
             log.Surv = log.Surv, Int = Int)
    }
    logPost.Bs.gammas2 <- function () {
        log.ptb <- event * drop(W2 %*% Bs.gammas) - log.Surv
        log.prior <- log.prior.Bs.gammas(Bs.gammas)
        sum(log.ptb, na.rm = TRUE) + log.prior
    }
    ArankDiff <- priorA.tauBs + 0.5 * qr(priorTau.Bs.gammas)$rank
    logPost.alphas <- function (alphas) {
        if (!paramRE) {
            Ms <- numeric(ns)
            if (estimateAlphas) {
                Mtime.alphas <- if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
            }
            if (paramExtra) {
                Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
            }
            log.Surv <- Int <- P * fastSumID(w * exp(log.h0s + Ms), id.GK)
            if (notNullW)
                log.Surv <- expWgammas * log.Surv
            log.ptb <- event * Mtime.alphas - log.Surv
            log.prior <- log.prior.alphas(alphas)
            list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv, Ms = Ms,
                 Int = Int)
        } else {
            Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
            log.Surv <- exp(Mtime) * Int
            if (notNullW)
                log.Surv <- expWgammas * log.Surv
            log.ptb <- event * Mtime - log.Surv
            log.prior <- log.prior.alphas(alphas)
            list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv)
        }
    }
    logPost.alphas2 <- function () {
        if (!paramRE) {
            Mtime.alphas <- if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        } else {
            Mtime.alphas <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
        }
        log.ptb <- event * Mtime.alphas - log.Surv
        log.prior <- log.prior.alphas(alphas)
        sum(log.ptb, na.rm = TRUE) + log.prior
    }
    logPost.Dalphas <- function (Dalphas) {
        Ms <- numeric(ns)
        if (estimateAlphas) {
            Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        }
        if (paramExtra) {
            Mtime.Dalphas <- if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
        }
        log.Surv <- Int <- P * fastSumID(w * exp(log.h0s + Ms), id.GK)
        if (notNullW)
            log.Surv <- expWgammas * log.Surv
        log.ptb <- event * Mtime.Dalphas - log.Surv
        log.prior <- log.prior.Dalphas(Dalphas)
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv,
             Ms = Ms, Int = Int)
    }
    logPost.Dalphas2 <- function () {
        Mtime.Dalphas <- if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
        log.ptb <- event * Mtime.Dalphas - log.Surv
        log.prior <- log.prior.Dalphas(Dalphas)
        sum(log.ptb, na.rm = TRUE) + log.prior
    }
    logPost.shape <- function (shape, which) {
        shapes[which] <- shape
        Ms <- numeric(ns)
        ###
        wFun <- w * weightFun(u.idGK, shapes, max.time)
        vl <- transFun.value(P * fastSumID(wFun * XsbetasZsb, id.GK), data.id)
        Mtime <- if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        ###
        wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
        vls <- transFun.value(P2 * fastSumID(wFun2 * XubetasZub, id.GK2), data.s)
        Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        ###
        if (paramExtra) {
            Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
        }
        ###
        log.Surv <- Int <- P * fastSumID(w * exp(log.h0s + Ms), id.GK)
        if (notNullW)
            log.Surv <- expWgammas * log.Surv
        log.ptb <- event * Mtime - log.Surv
        log.prior <- switch(which, "1" = log.prior.shape1(shape), 
                            "2" = log.prior.shape2(shape), "3" = log.prior.shape3(shape))
        list(log.post = sum(log.ptb, na.rm = TRUE) + log.prior, log.Surv = log.Surv,
             Ms = Ms, Int = Int, wFun = wFun, wFun2 = wFun2, vl = vl, vls = vls)
    }
    # define proposals
    # betas
    propCov.betas <- eigen(Covs$betas, symmetric = TRUE)
    scale.betas <- if (!is.null(ss <- scales[["betas"]])) ss else 5.66/nbetas
    r.betas <- function (n) {
        propCov.betas$values <- propCov.betas$values * scale.betas
        rmvnorm(n, mu = NULL, Sigma = propCov.betas)
    }
    # b
    propCov.RE <- lapply(Covs$b, eigen, symmetric = TRUE)
    scale.RE <- if (!is.null(ss <- scales[["b"]])) ss else rep(5.66/nRE, n)
    r.RE <- function (N) {
        out <- array(0, c(dim(b), N))
        for (i in 1:n) {
            propCov.RE[[i]][["values"]] <- propCov.RE[[i]][["values"]] * scale.RE[i]
            out[i, , ] <- rmvnorm(N, mu = NULL, Sigma = propCov.RE[[i]])
        }
        out
    }
    # invD
    isNulldf.RE <- is.null(df.RE)
    diagB <- diag(1, ncZ)
    if (isNulldf.RE) {
        K.invDn <- priorK.invD + n
        r.invD <- function (N) {
            drop(rWishart(N, K.invDn, R.Dbtb))
        }
    } else {
        K.invDn <- (df.RE / (df.RE - 2)) * (priorK.invD + n)
        r.invD <- function (N) {
            drop(rWishart(N, K.invDn, invD / K.invDn))
        }
    }
    # gammas
    if (notNullW) {
        propCov.gammas <- eigen(Covs$gammas, symmetric = TRUE)
        scale.gammas <- if (!is.null(ss <- scales$gammas)) ss else 5.66/ngammas
        r.gammas <- function (N) {
            propCov.gammas$values <- propCov.gammas$values * scale.gammas
            rmvnorm(N, mu = NULL, Sigma = propCov.gammas)
        }
    }
    # Bs.gammas
    propCov.Bs.gammas <- eigen(Covs$Bs.gammas, symmetric = TRUE)
    scale.Bs.gammas <- if (!is.null(ss <- scales$Bs.gammas)) ss else 5.66/nBs.gammas
    r.Bs.gammas <- function (N) {
        propCov.Bs.gammas$values <- propCov.Bs.gammas$values * scale.Bs.gammas
        rmvnorm(N, mu = NULL, Sigma = propCov.Bs.gammas)
    }
    # alphas
    if (estimateAlphas) {
        propCov.alphas <- eigen(Covs$alphas, symmetric = TRUE)
        scale.alphas <- if (!is.null(ss <- scales$alphas)) ss else 5.66/nalphas
        r.alphas <- function (N) {
            propCov.alphas$values <- propCov.alphas$values * scale.alphas
            rmvnorm(N, mu = NULL, Sigma = propCov.alphas)
        }
    }
    # Dalphas
    if (paramExtra) {
        propCov.Dalphas <- eigen(Covs$Dalphas, symmetric = TRUE)
        scale.Dalphas <- if (!is.null(ss <- scales$Dalphas)) ss else 5.66/nDalphas
        r.Dalphas <- function (N) {
            propCov.Dalphas$values <- propCov.Dalphas$values * scale.Dalphas
            rmvnorm(N, mu = NULL, Sigma = propCov.Dalphas)
        }
    }
    # number of iterations
    n.adapt <- control$n.adapt
    n.burnin <- control$n.burnin
    totalIter <- control$n.iter + n.adapt + n.burnin
    n.thin <- control$n.thin
    n.batch <- control$n.batch
    # objects to keep results
    resInd <- seq(n.adapt + n.burnin + 1L, totalIter, by = n.thin)
    n.out <- length(resInd)
    res.betas <- matrix(0, n.out, length(betas))
    if (hasScale)
        res.tau <- matrix(0, n.out, 1)
    res.b <- array(0, c(dim(b), n.out))
    res.invD <- matrix(0, n.out, length(invD))
    res.Bs.gammas <- matrix(0, n.out, length(Bs.gammas))
    if (baseHazP)
        res.tauBs <- matrix(0, n.out, 1)
    if (notNullW)
        res.gammas <- matrix(0, n.out, length(gammas))
    if (estimateAlphas)
        res.alphas <- matrix(0, n.out, length(alphas))
    if (paramExtra)
        res.Dalphas <- matrix(0, n.out, length(Dalphas))
    if (estimateWeightFun)
        res.shapes <- matrix(0, n.out, length(shapes))
    res.logLik <- matrix(0, n.out, n)
    # acceptance rates
    ar.betas <- ar.invD <- ar.gammas <- ar.Bs.gammas <- ar.alphas <- ar.Dalphas <- numeric(totalIter)
    ar.b <- matrix(0, totalIter, n)
    # initiate all components at the starting values
    Xbetas <- drop(X %*% betas)
    Zb <- rowSums(Z * b[id, , drop = FALSE])
    eta.y <- Xbetas + Zb
    Mtime <- numeric(n)
    Ms <- numeric(ns)
    if (paramValue) {
        Xtimebetas <- drop(Xtime %*% betas)
        Ztimeb <- rowSums(Ztime * b)
        Xsbetas <- drop(Xs %*% betas)
        Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
        vl <- transFun.value(Xtimebetas + Ztimeb, data.id)
        vls <- transFun.value(Xsbetas + Zsb, data.s)
        is.matrix.vl <- is.matrix(vl); is.matrix.vls <- is.matrix(vls)
        Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
    }
    if (paramExtra) {
        Xtime.extrabetas <- drop(Xtime.extra %*% betas[iF])
        Ztime.extrab <- rowSums(Ztime.extra * b[, iR, drop = FALSE])
        Xs.extrabetas <- drop(Xs.extra %*% betas[iF])
        Zs.extrab <- rowSums(Zs.extra * b[id.GK, iR, drop = FALSE])
        ex <- transFun.extra(Xtime.extrabetas + Ztime.extrab, data.id)
        exs <- transFun.extra(Xs.extrabetas + Zs.extrab, data.s)
        is.matrix.ex <- is.matrix(ex); is.matrix.exs <- is.matrix(exs)
        Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
        Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
    }
    if (estimateWeightFun) {
        wFun <- w * weightFun(u.idGK, shapes, max.time)
        Xsbetas <- drop(Xs %*% betas)
        Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
        wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
        Xubetas <- drop(Xu %*% betas)
        Zub <- rowSums(Zu * b[id.GKu, , drop = FALSE])
        vl <- transFun.value(P * fastSumID(wFun * (Xsbetas + Zsb), id.GK), data.id)
        vls <- transFun.value(P2 * fastSumID(wFun2 * (Xubetas + Zub), id.GK2), data.s)
        is.matrix.vl <- is.matrix(vl); is.matrix.vls <- is.matrix(vls)
        Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
        Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
    }
    log.h0s <- drop(W2s %*% Bs.gammas)
    log.Surv <- Int <- P * if (!paramRE) {
        fastSumID(w * exp(log.h0s + Ms), id.GK)
    } else {
        fastSumID(w * exp(log.h0s), id.GK)
    }
    if (paramRE) {
        Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
        log.Surv <- log.Surv * exp(Mtime)
    }
    if (notNullW) {
        expWgammas <- exp(drop(W %*% gammas))
        log.Surv <- log.Surv * expWgammas
    }
    #########################################################################################################
    # run the MCMC
    set.seed(control$seed)
    batch <- 1L
    jj <- 0L
    batchStart <- round((n.adapt + n.burnin) / n.batch)
    if (control$verbose) {
        cat("\n MCMC iterations:\n\n")
        pb <- txtProgressBar(0, totalIter, style = 3, char = "+", width = 50)
    }
    time <- system.time(for (i in seq_len(totalIter)) {
        if (i == 1L || !i %% n.batch) {
            if (i > 1 && i <= n.adapt) {
                ss <- seq(1L + n.batch * (batch - 1L), n.batch * batch)
                if (is.null(scales[["betas"]]))
                    scale.betas <- scbetasF <- adjustScaleRW(scale.betas, mean(ar.betas[ss]), nbetas)
                if (is.null(scales[["b"]]))
                    scale.RE <- mapply(adjustScaleRW, scale.RE, colMeans(ar.b[ss, ]), nRE)
                if (notNullW && is.null(scales$gammas)) {
                    scale.gammas <- scgammasF <- adjustScaleRW(scale.gammas, mean(ar.gammas[ss]), ngammas)
                }
                if (is.null(scales$Bs.gammas))
                    scale.Bs.gammas <- scBs.gammasF <- adjustScaleRW(scale.Bs.gammas, mean(ar.Bs.gammas[ss]), nBs.gammas)
                if (estimateAlphas && is.null(scales$alphas)) {
                    scale.alphas <- scalphasF <- adjustScaleRW(scale.alphas, mean(ar.alphas[ss]), nalphas)
                }
                if (paramExtra && is.null(scales$Dalphas)) {
                    scale.Dalphas <- adjustScaleRW(scale.Dalphas, mean(ar.Dalphas[ss]), nDalphas)
                }
                if (!isNulldf.RE) {
                    K.invDn <- adjustKRW(K.invDn, mean(ar.invD[ss]), ncZ)
                }
            }
            if (i > 1)
                batch <- batch + 1L
            new.betas <- r.betas(n.batch)
            new.b <- r.RE(n.batch)
            if (notNullW)
                new.gammas <- r.gammas(n.batch)
            new.Bs.gammas <- r.Bs.gammas(n.batch)
            if (estimateAlphas)
                new.alphas <- r.alphas(n.batch)
            if (paramExtra)
                new.Dalphas <- r.Dalphas(n.batch)
        }
        # batch index
        ii <- i - n.batch * if (!i %% n.batch) i %/% n.batch - 1L else i %/% n.batch
        # update betas
        lP.old.betas <- logPost.betas2()
        new.betas[ii, ] <- new.betas[ii, ] + betas
        lP.betas <- logPost.betas(new.betas[ii, ])
        lP.new.betas <- lP.betas$log.post
        lRatio.betas <- lP.new.betas - lP.old.betas
        if (lRatio.betas >= 0 || runif(1L) < exp(lRatio.betas)) {
            ar.betas[i] <- 1
            betas <- new.betas[ii, ]
            Xbetas <- lP.betas$Xbetas
            eta.y <- Xbetas + Zb
            if (!paramRE) {
                Xtimebetas <- lP.betas$Xtimebetas
                Xsbetas <- lP.betas$Xsbetas
                Xtime.extrabetas <- lP.betas$Xtime.extrabetas
                Xs.extrabetas <- lP.betas$Xs.extrabetas
                Xubetas <- lP.betas$Xubetas
                vl <- lP.betas$vl; vls <- lP.betas$vls
                ex <- lP.betas$ex; exs <- lP.betas$exs
                Mtime <- lP.betas$Mtime
                Ms <- lP.betas$Ms
                log.Surv <- lP.betas$log.Surv
                Int <- lP.betas$Int
            }
            if (param == "shared-betasRE")
                log.Surv <- lP.betas$log.Surv
        }
        # update tau
        if (hasScale) {
            tau <- slice.tau(logPost.tau, tau, step = 0.5)
        }
        # update RE
        lP.old.b <- logPost.RE2()
        new.b[, , ii] <- new.b[, , ii] + b
        lP.RE <- logPost.RE(as.matrix(new.b[, , ii]))
        lP.new.b <- lP.RE$log.post
        lRatio.b <- lP.new.b - lP.old.b
        indRE <- runif(n) < pmin(exp(lRatio.b), 1)
        indRE.GK <- indRE[id.GK]
        indRE.id <- indRE[id]
        ar.b[i, indRE] <- 1
        b[indRE, ] <- new.b[indRE, , ii]
        Zb[indRE.id] <- lP.RE$Zb[indRE.id]
        log.Surv[indRE] <- lP.RE$log.Surv[indRE]
        eta.y <- Xbetas + Zb
        if (!paramRE) {
            Int[indRE] <- lP.RE$Int[indRE]
            Ms[indRE.GK] <- lP.RE$Ms[indRE.GK]
            if (estimateAlphas) {
                if (notestimateWeightFun)
                    Ztimeb[indRE] <- lP.RE$Ztimeb[indRE]
                Zsb[indRE.GK] <- lP.RE$Zsb[indRE.GK]
                if (is.matrix.vl) {
                    vl[indRE, ] <- lP.RE$vl[indRE, ]
                    vls[indRE.GK, ] <- lP.RE$vls[indRE.GK, ]
                } else {
                    vl[indRE] <- lP.RE$vl[indRE]
                    vls[indRE.GK] <- lP.RE$vls[indRE.GK]
                }
            }
            if (paramExtra) {
                Ztime.extrab[indRE] <- lP.RE$Ztime.extrab[indRE]
                Zs.extrab[indRE.GK] <- lP.RE$Zs.extrab[indRE.GK]
                if (is.matrix.ex) {
                    ex[indRE, ] <- lP.RE$ex[indRE, ]
                    exs[indRE.GK, ] <- lP.RE$exs[indRE.GK, ]
                } else {
                    ex[indRE] <- lP.RE$ex[indRE]
                    exs[indRE.GK] <- lP.RE$exs[indRE.GK]
                }
            }
            if (estimateWeightFun) {
                indRE.GKu <- indRE[id.GKu]
                Zub[indRE.GKu] <- lP.RE$Zub[indRE.GKu]
            }
        }
        # update invD
        if (isNulldf.RE) {
            R.Dbtb <- solve.default(priorR.invD + crossprod(b), diagB)
            new.invD <- r.invD(1)
            ar.invD[i] <- 1
            invD <- new.invD
        } else {
            new.invD <- r.invD(1)
            lP.old.invD <- logPost.invD(invD)
            lP.new.invD <- logPost.invD(new.invD)
            lRatio.invD <- lP.new.invD + dwish(invD, invD / K.invDn, K.invDn, TRUE) -
                lP.old.invD - dwish(new.invD, invD / K.invDn, K.invDn, TRUE)
            if (lRatio.invD >= 0 || runif(1L) < exp(lRatio.invD)) {
                ar.invD[i] <- 1
                invD <- new.invD
            }
        }
        # update gammas
        if (notNullW) {
            lP.old.gammas <- logPost.gammas(gammas)$log.post
            new.gammas[ii, ] <- new.gammas[ii, ] + gammas
            lP.gammas <- logPost.gammas(new.gammas[ii, ])
            lP.new.gammas <- lP.gammas$log.post
            lRatio.gammas <- lP.new.gammas - lP.old.gammas
            if (lRatio.gammas >= 0 || runif(1L) < exp(lRatio.gammas)) {
                ar.gammas[i] <- 1
                gammas <- new.gammas[ii, ]
                expWgammas <- lP.gammas$expWgammas
                log.Surv <- lP.gammas$log.Surv
            }
        }
        # update Bs.gammas
        lP.old.Bs.gammas <- logPost.Bs.gammas2()
        new.Bs.gammas[ii, ] <- new.Bs.gammas[ii, ] + Bs.gammas
        lP.Bs.gammas <- logPost.Bs.gammas(new.Bs.gammas[ii, ])
        lP.new.Bs.gammas <- lP.Bs.gammas$log.post
        lRatio.Bs.gammas <- lP.new.Bs.gammas - lP.old.Bs.gammas
        if (lRatio.Bs.gammas >= 0 || runif(1L) < exp(lRatio.Bs.gammas)) {
            ar.Bs.gammas[i] <- 1
            Bs.gammas <- new.Bs.gammas[ii, ]
            log.h0s <- lP.Bs.gammas$log.h0s
            log.Surv <- lP.Bs.gammas$log.Surv
            Int <- lP.Bs.gammas$Int
        }
        # update tauBs
        if (baseHazP) {
            BB <- priorB.tauBs + 0.5 * drop(crossprod(Bs.gammas, priorTau.Bs.gammas %*% Bs.gammas))
            tauBs <- rgamma(1L, ArankDiff, BB)
        }
        # update alphas
        if (estimateAlphas) {
            lP.old.alphas <- logPost.alphas2()
            new.alphas[ii, ] <- new.alphas[ii, ] + alphas
            lP.alphas <- logPost.alphas(new.alphas[ii, ])
            lP.new.alphas <- lP.alphas$log.post
            lRatio.alphas <- lP.new.alphas - lP.old.alphas
            if (lRatio.alphas >= 0 || runif(1L) < exp(lRatio.alphas)) {
                ar.alphas[i] <- 1
                alphas <- new.alphas[ii, ]
                log.Surv <- lP.alphas$log.Surv
                if (!paramRE) {
                    Ms <- lP.alphas$Ms
                    Int <- lP.alphas$Int
                }
            }
        }
        # update Dalphas
        if (paramExtra) {
            lP.old.Dalphas <- logPost.Dalphas2()
            new.Dalphas[ii, ] <- new.Dalphas[ii, ] + Dalphas
            lP.Dalphas <- logPost.Dalphas(new.Dalphas[ii, ])
            lP.new.Dalphas <- lP.Dalphas$log.post
            lRatio.Dalphas <- lP.new.Dalphas - lP.old.Dalphas
            if (lRatio.Dalphas >= 0 || runif(1L) < exp(lRatio.Dalphas)) {
                ar.Dalphas[i] <- 1
                Dalphas <- new.Dalphas[ii, ]
                Ms <- lP.Dalphas$Ms
                log.Surv <- lP.Dalphas$log.Surv
                Int <- lP.Dalphas$Int
            }
        }
        # update shapes
        if (estimateWeightFun) {
            if (control$verbose2)
                cat("\ni =", i, "\tshapes =", round(shapes, 3), 
                    if (paramExtra) "\tDalphas =", if (paramExtra) round(Dalphas, 3), 
                    "\talphas =", round(alphas, 3), "\tbetas = ", round(betas, 3))
            XsbetasZsb <- Xsbetas + Zsb
            XubetasZub <- Xubetas + Zub
            for (shp in seq.nshapes) {
                ss <- 1
                slice.shp <- slice.shape(logPost.shape, shapes, step = ss, which = shp)
                while(slice.shp$fail) {
                    ss <- ss/10
                    if (ss < 1e-05)
                        break
                    slice.shp <- slice.shape(logPost.shape, shapes, step = ss, which = shp)
                }
                shapes[shp] <- slice.shp$new.shape
                if (shp == nshapes) {
                    log.Surv <- slice.shp$log.Surv
                    Ms <- slice.shp$Ms
                    Int <- slice.shp$Int
                    wFun <- slice.shp$wFun
                    wFun2 <- slice.shp$wFun2
                    vl <- slice.shp$vl
                    vls <- slice.shp$vls                    
                }
            }
        }
        if (control$verbose && !i %% n.batch)
            setTxtProgressBar(pb, i)
        # adapt
        if (control$adapt && jj > 300 && !i %% n.batch) {
            ss <- seq(1 + n.batch * (batch - 2), n.batch * (batch - 1))
            Covs$betas <- 0.8*Covs$betas + 0.2*var(res.betas)
            propCov.betas <- eigen(Covs$betas, symmetric = TRUE)
            scale.betas <- adjustScaleRW(scale.betas, mean(ar.betas[ss]), nbetas, batch - batchStart, scbetasF)
            Covs$Bs.gammas <- 0.8*Covs$Bs.gammas + 0.2*var(res.Bs.gammas)
            propCov.Bs.gammas <- eigen(Covs$Bs.gammas, symmetric = TRUE)
            scale.Bs.gammas <- adjustScaleRW(scale.Bs.gammas, mean(ar.Bs.gammas[ss]), nBs.gammas,
                                             batch - batchStart, scBs.gammasF)
            if (notNullW) {
                Covs$gammas <- 0.8*Covs$gammas + 0.2*var(res.gammas)
                propCov.gammas <- eigen(Covs$gammas, symmetric = TRUE)
                scale.gammas <- adjustScaleRW(scale.gammas, mean(ar.gammas[ss]), ngammas, batch - batchStart, scgammasF)
            }
            if (estimateAlphas) {
                Covs$alphas <- 0.8*Covs$alphas + 0.2*var(res.alphas)
                propCov.alphas <- eigen(Covs$alphas, symmetric = TRUE)
                scale.alphas <- adjustScaleRW(scale.alphas, mean(ar.alphas[ss]), nalphas, batch - batchStart, scalphasF)
            }
        }
        # save results
        if (i %in% resInd) {
            jj <- match(i, resInd)
            res.betas[jj, ] <- betas
            if (hasScale)
                res.tau[jj, ] <- tau
            res.b[, , jj] <- b
            res.invD[jj, ] <- c(invD)
            if (notNullW)
                res.gammas[jj, ] <- gammas
            res.Bs.gammas[jj, ] <- Bs.gammas
            if (baseHazP)
                res.tauBs[jj, ] <- tauBs
            if (estimateAlphas)
                res.alphas[jj, ] <- alphas
            if (paramExtra)
                res.Dalphas[jj, ] <- Dalphas
            if (estimateWeightFun)
                res.shapes[jj, ] <- shapes
            log.pyb <- fastSumID(densLong(y.long, eta.y, 1/sqrt(tau), log = TRUE, data), id)
            log.h0s <- drop(W2s %*% Bs.gammas)
            if (!paramRE) {
                Mtime <- numeric(n)
                Ms <- numeric(ns)
                if (paramValue) {
                    Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                    Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
                }
                if (paramExtra) {
                    Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
                    Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
                }
                if (estimateWeightFun) {
                    wFun <- w * weightFun(u.idGK, shapes, max.time)
                    Xsbetas <- drop(Xs %*% betas)
                    Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
                    vl <- transFun.value(P * fastSumID(wFun * (Xsbetas + Zsb), id.GK), data.id)
                    Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
                    wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
                    Xubetas <- drop(Xu %*% betas)
                    Zub <- rowSums(Zu * b[id.GKu, , drop = FALSE])
                    vls <- transFun.value(P2 * fastSumID(wFun2 * (Xubetas + Zub), id.GK2), data.s)
                    Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
                }
                log.Surv <- Int <- P * fastSumID(w * exp(log.h0s + Ms), id.GK)
                if (notNullW)
                    log.Surv <- expWgammas * log.Surv
            } else {
                log.Surv <- Int <- P * fastSumID(w * exp(log.h0s), id.GK)
                Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
                log.Surv <- exp(Mtime) * log.Surv
                if (notNullW)
                    log.Surv <- expWgammas * log.Surv
            }
            log.h <- drop(W2 %*% Bs.gammas) + Mtime
            if (notNullW)
                log.h <- drop(W %*% gammas) + log.h
            log.ptb <- event * log.h - log.Surv
            log.pb <- densRE(b, invD = invD, log = TRUE, prop = FALSE)
            res.logLik[jj, ] <- log.pyb + log.ptb + log.pb
        }
    })
    if (control$verbose)
        close(pb)    
    mcmcOut <- list(betas = res.betas, sigma = if (hasScale) 1/sqrt(res.tau), b = res.b,
                    D = if (ncZ > 1) t(apply(res.invD, 1, function (x) solve.default(matrix(x, ncZ))))
                    else as.matrix(apply(res.invD, 1, function (x) solve.default(matrix(x, ncZ)))),
                    gammas = if (notNullW) res.gammas, Bs.gammas = res.Bs.gammas,
                    tauBs = if (baseHazP) res.tauBs,
                    alphas = if (estimateAlphas) res.alphas, Dalphas = if (paramExtra) res.Dalphas,
                    shapes = if (estimateWeightFun) res.shapes)
    mcmcOut <- mcmcOut[!sapply(mcmcOut, is.null)]
    # calculate pD
    D.bar <- - 2 * mean(rowSums(res.logLik, na.rm = TRUE), na.rm = TRUE)
    postMeans <- lapply(mcmcOut, function (x) {
        d <- dim(x)
        if (!is.null(d) && length(d) > 2) apply(x, c(1, 2), mean) else colMeans(as.matrix(x))
    })
    dim(postMeans$D) <- c(ncZ, ncZ)
    betas <- postMeans$betas; sigma <- postMeans$sigma; b <- postMeans$b; D <- postMeans$D
    gammas <- postMeans$gammas; Bs.gammas <- postMeans$Bs.gammas; alphas <- postMeans$alphas;
    Dalphas <- postMeans$Dalphas; shapes <- postMeans$shapes
    log.pyb <- fastSumID(densLong(y.long, drop(X %*% betas) + rowSums(Z * b[id, , drop = FALSE]),
                                  sigma, log = TRUE, data), id)
    log.h0s <- drop(W2s %*% Bs.gammas)
    if (!paramRE) {
        Mtime <- numeric(n)
        Ms <- numeric(ns)
        if (paramValue) {
            vl <- transFun.value(drop(Xtime %*% betas) + rowSums(Ztime * b), data.id)
            vls <- transFun.value(drop(Xs %*% betas) + rowSums(Zs * b[id.GK, , drop = FALSE]), data.s)
            Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        }
        if (paramExtra) {
            ex <- transFun.extra(drop(Xtime.extra %*% betas[iF]) +
                                     rowSums(Ztime.extra * b[, iR, drop = FALSE]), data.id)
            exs <- transFun.extra(drop(Xs.extra %*% betas[iF]) +
                                      rowSums(Zs.extra * b[id.GK, iR, drop = FALSE]), data.s)
            Mtime <- Mtime + if (is.matrix.ex) drop(ex %*% Dalphas) else ex * Dalphas
            Ms <- Ms + if (is.matrix.exs) drop(exs %*% Dalphas) else exs * Dalphas
        }
        if (estimateWeightFun) {
            wFun <- w * weightFun(u.idGK, shapes, max.time)
            Xsbetas <- drop(Xs %*% betas)
            Zsb <- rowSums(Zs * b[id.GK, , drop = FALSE])
            vl <- transFun.value(P * fastSumID(wFun * (Xsbetas + Zsb), id.GK), data.id)
            Mtime <- Mtime + if (is.matrix.vl) drop(vl %*% alphas) else vl * alphas
            wFun2 <- w2 * weightFun(u.idGK2, shapes, max.time)
            Xubetas <- drop(Xu %*% betas)
            Zub <- rowSums(Zu * b[id.GKu, , drop = FALSE])
            vls <- transFun.value(P2 * fastSumID(wFun2 * (Xubetas + Zub), id.GK2), data.s)
            Ms <- Ms + if (is.matrix.vls) drop(vls %*% alphas) else vls * alphas
        }
        log.Surv <- P * fastSumID(w * exp(log.h0s + Ms), id.GK)
        if (notNullW)
            log.Surv <- exp(drop(W %*% gammas)) * log.Surv
    } else {
        Mtime <- if (paramSharedRE) drop(b %*% alphas) else drop((betas[indBetas2] + b) %*% alphas)
        log.Surv <- exp(Mtime) * P * fastSumID(w * exp(log.h0s), id.GK)
        if (notNullW)
            log.Surv <- exp(drop(W %*% gammas)) * log.Surv
    }
    log.h <- drop(W2 %*% Bs.gammas) + Mtime
    if (notNullW)
        log.h <- drop(W %*% gammas) + log.h
    log.ptb <- event * log.h - log.Surv
    log.pb <- densRE(b, D = D, log = TRUE, prop = FALSE)
    D.hat <- - 2 * sum(log.pyb + log.ptb + log.pb, na.rm = TRUE)
    pD <- D.bar - D.hat
    indb <- names(mcmcOut) != "b"
    postVarsRE <- apply(res.b, 1, function (x) var(t(x)))
    dim(postVarsRE)<- c(ncZ, ncZ, n)
    keepD <- length(betas) + 1 + which(!lower.tri(invD, TRUE))
    keepAR <- -seq_len(n.adapt)
    postModes <- lapply(mcmcOut[indb], function (x) apply(as.matrix(x), 2, modes))
    dim(postModes$D) <- c(ncZ, ncZ)
    list(mcmc = if (control$keepRE) mcmcOut else mcmcOut[indb], postMeans = postMeans,
         postModes = postModes,
         postVarsRE = postVarsRE,
         StErr = lapply(mcmcOut[indb], stdErr),
         EffectiveSize = lapply(mcmcOut[indb], effectiveSize),
         StDev = lapply(mcmcOut[indb], function (x) apply(as.matrix(x), 2, sd)),
         CIs = lapply(mcmcOut[indb], function (x) apply(as.matrix(x), 2, quantile, probs = c(0.025, 0.975))),
         Pvalues = lapply(mcmcOut[indb], function (x) apply(as.matrix(x), 2, computeP)),
         vcov = if (ncZ > 1) var(do.call(cbind, mcmcOut[indb])[, -keepD]) else var(do.call(cbind, mcmcOut[indb])),
         pD = pD, DIC = pD + D.bar, CPO = 1 / colMeans(exp(-res.logLik)),
         LPML = sum(log(1 / colMeans(exp(-res.logLik))), na.rm = TRUE), time = time,
         scales = list(betas = scale.betas, b = scale.RE, Bs.gammas = scale.Bs.gammas,
                       gammas = if (notNullW) scale.gammas,
                       alphas = if (estimateAlphas) scale.alphas,
                       Dalphas = if (paramExtra) scale.Dalphas),
         Covs = Covs,
         acceptRates = list(betas = mean(ar.betas[keepAR]), b = colMeans(ar.b[keepAR, ]),
                            D = mean(ar.invD[keepAR]), Bs.gammas = mean(ar.Bs.gammas[keepAR]),
                            gammas = if (notNullW) mean(ar.gammas[keepAR]),
                            alphas = if (estimateAlphas) mean(ar.alphas[keepAR]),
                            Dalphas = if (paramExtra) mean(ar.Dalphas[keepAR])))
}
