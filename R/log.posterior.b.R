log.posterior.b <-
function (b, y, Mats, ii) {
    id.i <- id %in% ii
    idT.i <- idT %in% ii
    ids.i <- ids %in% ii
    X.i <- X[id.i, , drop = FALSE]
    Z.i <- Z[id.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% betas.new + Z.i %*% b)
    logY <- densLong(y[id.i], mu.y, sigma.new, log = TRUE, data = newdata[id.i, ])
    log.p.yb <- sum(logY)
    log.p.b <- densRE(b, D = D.new, log = TRUE, prop = FALSE)
    st <- Mats[[ii]]$st
    wk <- Mats[[ii]]$wk
    P <- Mats[[ii]]$P
    W2s <- Mats[[ii]]$W2s
    Xs <- Mats[[ii]]$Xs
    Zs <- Mats[[ii]]$Zs
    Xs.extra <- Mats[[ii]]$Xs.extra
    Zs.extra <- Mats[[ii]]$Zs.extra
    ind <- Mats[[ii]]$ind
    idT <- Mats[[ii]]$idT
    if (param %in% c("td-value", "td-both"))
        Ys <- transFun.value(c(Xs %*% betas.new + Zs %*% b), data.s[ids.i, ])
    if (param %in% c("td-extra", "td-both"))
        Ys.extra <- transFun.extra(c(Xs.extra %*% betas.new[indFixed] + Zs.extra %*% b[indRandom]), data.s[ids.i, ])
    tt <- c(switch(param,
                   "td-value" = as.matrix(Ys) %*% alphas.new, 
                   "td-extra" =  as.matrix(Ys.extra) %*% Dalphas.new,
                   "td-both" = as.matrix(Ys) %*% alphas.new + as.matrix(Ys.extra) %*% Dalphas.new,
                   "shared-betasRE" = rep(sum((betas[indBetas] + b) * alphas.new), length(st)),
                   "shared-RE" = rep(sum(b * alphas.new), length(st))))
    eta.tw <- if (!is.null(W)) {
            as.vector(W[ii, , drop = FALSE] %*% gammas.new)
    } else 0
    Vi <- exp(c(W2s %*% Bs.gammas.new) + tt)
    log.survival <- - sum(exp(eta.tw) * P * fastSumID(wk * Vi, idT))
    if (all(st == 0))
        log.survival <- 1
    log.p.yb + log.survival + log.p.b
}
