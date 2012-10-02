S.b <-
function (t, b, ii, Mats) {
    if (t == 0)
        return(1)
    idT.i <- idT %in% ii
    st <- Mats$st
    wk <- Mats$wk
    P <- Mats$P
    Xs <- Mats$Xs
    Zs <- Mats$Zs
    Xs.deriv <- Mats$Xs.deriv
    Zs.deriv <- Mats$Zs.deriv
    Ws.intF.vl <- Mats$Ws.intF.vl
    Ws.intF.sl <- Mats$Ws.intF.sl
    ind <- Mats$ind
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
    exp(log.survival)
}
