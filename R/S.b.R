S.b <-
function (t, b, ii, Mats) {
    if (t == 0)
        return(1)
    idT.i <- idT %in% ii
    ids.i <- ids %in% ii
    st <- Mats$st
    wk <- Mats$wk
    P <- Mats$P
    Xs <- Mats$Xs
    Zs <- Mats$Zs
    Xs.extra <- Mats$Xs.extra
    Zs.extra <- Mats$Zs.extra
    W2s <- Mats$W2s
    ind <- Mats$ind
    idT <- Mats$idT
    if (param %in% c("td-value", "td-both"))
        Ys <- transFun.value(c(Xs %*% betas.new + Zs %*% b), data.s[ids.i, ])
    if (param %in% c("td-extra", "td-both"))
        Ys.extra <- transFun.extra(c(Xs.extra %*% betas.new[indFixed] + 
                                         Zs.extra %*% b[indRandom]), data.s[ids.i, ])
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
    ind <- Mats$st < min(Mats$kn)
    wk[ind] <- 0
    log.survival <- - sum(exp(eta.tw) * P * fastSumID(wk * Vi, idT))
    exp(log.survival)
}
