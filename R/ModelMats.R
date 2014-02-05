ModelMats <-
function (time, ii) {
    id.GK <- rep(ii, each = 15)
    GQsurv <- if (object$control$GQsurv == "GaussKronrod") gaussKronrod() else gaussLegendre(object$control$GQsurv.k)
    wk <- GQsurv$wk
    sk <- GQsurv$sk
    P <- time / 2
    st <- P * (sk + 1)
    data.id2 <- data.id[id.GK, ]
    data.id2[[timeVar]] <- pmax(st - lag, 0)
    kn <- object$control$knots
    W2s <- splineDesign(unlist(kn, use.names = FALSE), st, 
                        ord = object$control$ordSpline, outer.ok = TRUE)    
    out <- list(st = st, wk = rep(wk, length(P)), P = P, W2s = W2s, kn = kn, 
                idT = rep(seq_along(P), each = 15))
    if (param %in% c("td-value", "td-both")) {
        mfX <- model.frame(delete.response(TermsX), data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        out$Xs <- model.matrix(formYx, mfX)
        out$Zs <- model.matrix(formYz, mfZ)
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.extra <- model.frame(TermsX.extra, data = data.id2)
        mfZ.extra <- model.frame(TermsZ.extra, data = data.id2)
        out$Xs.extra <- model.matrix(extraForm$fixed, mfX.extra)
        out$Zs.extra <- model.matrix(extraForm$random, mfZ.extra)
    }
    out
}
