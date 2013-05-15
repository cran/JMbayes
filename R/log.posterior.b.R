log.posterior.b <-
function (b, y, Mats, survMod, ii) {
    id.i <- id %in% ii
    idT.i <- idT %in% ii
    X.i <- X[id.i, , drop = FALSE]
    Z.i <- Z[id.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% betas.new) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
    logY <- if (!robust) dnorm(y[id.i], mu.y, sigma.new, TRUE) else dgt(y[id.i], mu.y, sigma.new, df, TRUE)
    log.p.yb <- sum(logY)
    log.p.b <- if (!robust.b) dmvnorm(b, rep(0, ncol(Z)), D.new, TRUE) else dmvt(b, rep(0, ncol(Z)), D.new, df.b, TRUE)
    st <- Mats[[ii]]$st
    wk <- Mats[[ii]]$wk
    P <- Mats[[ii]]$P
    Xs <- Mats[[ii]]$Xs
    Zs <- Mats[[ii]]$Zs
    Xs.deriv <- Mats[[ii]]$Xs.deriv
    Zs.deriv <- Mats[[ii]]$Zs.deriv
    Ws.intF.vl <- Mats[[ii]]$Ws.intF.vl
    Ws.intF.sl <- Mats[[ii]]$Ws.intF.sl
    ind <- Mats[[ii]]$ind
    if (param %in% c("td-value", "td-both"))
        Ys <- as.vector(Xs %*% betas.new + rowSums(Zs * rep(b, each = nrow(Zs))))
    if (param %in% c("td-extra", "td-both"))
        Ys.deriv <- as.vector(Xs.deriv %*% betas.new[indFixed]) + 
            rowSums(Zs.deriv * rep(b[indRandom], each = nrow(Zs)))
    tt <- switch(param,
        "td-value" = c(Ws.intF.vl %*% alpha.new) * Ys, 
        "td-extra" = c(Ws.intF.sl %*% Dalpha.new) * Ys.deriv,
        "td-both" = c(Ws.intF.vl %*% alpha.new) * Ys + 
            c(Ws.intF.sl %*% Dalpha.new) * Ys.deriv,
        "shared-RE" = rep(sum(b * alpha.new), length(st)))
    eta.tw <- if (!is.null(W)) {
            as.vector(W[ii, , drop = FALSE] %*% gammas.new)
    } else 0
    log.survival <- if (survMod == "weibull-PH") {
        Vi <- exp(log(sigma.t.new) + (sigma.t.new - 1) * log(st) + tt)
        - exp(eta.tw) * P * sum(wk * Vi)
    } else if (survMod == "spline-PH") {
        kn <- object$control$knots
        W2s <- splineDesign(unlist(kn, use.names = FALSE), st, 
            ord = object$control$ordSpline, outer.ok = TRUE)
        Vi <- exp(c(W2s %*% Bs.gammas.new) + tt)
        idT <- rep(seq_along(P), each = 15)
        - sum(exp(eta.tw) * P * tapply(wk * Vi, idT, sum))
    }
    if (all(st == 0))
        log.survival <- 1
    log.p.yb + log.survival + log.p.b
}
