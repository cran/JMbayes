fitted.JMbayes <-
function (object, process = c("Longitudinal", "longitudinal", "Event", "event"), 
                            type = c("Marginal", "marginal", "Subject", "subject"), nullY = FALSE, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    process <- match.arg(process)
    type <- match.arg(type)
    if (process == "Longitudinal" || process == "longitudinal") {
        fitY <- c(object$x$X %*% object$coefficients$betas)
        names(fitY) <- names(object$y$y)
        if (type == "Subject" || type == "subject")
            fitY <- fitY + rowSums(object$x$Z * ranef(object)[object$id, ])
        fitY
    } else {
        survMod <- object$survMod
        timeVar <- object$timeVar
        robust <- object$robust
        robust.b <- object$robust.b
        df <- object$df
        df.b <- object$df.b
        param <- object$param
        extraForm <- object$extraForm
        indFixed <- extraForm$indFixed
        indRandom <- extraForm$indRandom
        lag <- object$y$lag
        TermsX <- object$termsYx
        TermsZ <- object$termsYz
        TermsX.deriv <- object$termsYx.deriv
        TermsZ.deriv <- object$termsYz.deriv
        formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
        formYz <- object$formYz
        times <- object$data[[timeVar]]
        wk <- gaussKronrod()$wk
        sk <- gaussKronrod()$sk
        K <- length(sk)
        P <- times/2
        st <- outer(P, sk + 1)
        id.GK <- rep(seq_along(times), each = K)
        data.id2 <- object$data.id[rep(object$id, each = K), ]
        data.id2[[timeVar]] <- pmax(c(t(st)) - lag, 0)
        if (param %in% c("td-value", "td-both")) {
            mfX <- model.frame(TermsX, data = data.id2)
            mfZ <- model.frame(TermsZ, data = data.id2)
            Xs <- model.matrix(formYx, mfX)
            Zs <- model.matrix(formYz, mfZ)
        }
        if (param %in% c("td-extra", "td-both")) {
            mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
            mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
            Xs.deriv <- model.matrix(extraForm$fixed, mfX.deriv)
            Zs.deriv <- model.matrix(extraForm$random, mfZ.deriv)
        }
        betas <- object$coefficients$betas
        sigma <- object$coefficients$sigma
        D <- object$coefficients$D
        gammas <- object$coefficients$gammas
        alpha <- object$coefficients$alpha
        Dalpha <- object$coefficients$Dalpha
        if (nullY) {
            alpha <- rep(0, length.out = length(alpha))
            Dalpha <- 0
        }
        sigma.t <- object$coefficients$sigma.t
        Bs.gammas <- object$coefficients$Bs.gammas
        b <- ranef(object)
        idK <- rep(object$id, each = K)
        b <- b[idK, ]
        if (param %in% c("td-value", "td-both")) {
            Ys <- as.vector(Xs %*% betas + rowSums(Zs * b))
        }
        if (param %in% c("td-extra", "td-both")) {
            Ys.deriv <- as.vector(Xs.deriv %*% betas[indFixed]) + 
                rowSums(Zs.deriv * b[, indRandom, drop = FALSE])
        }
        tt <- switch(param,
                     "td-value" = alpha * Ys, 
                     "td-extra" = Dalpha * Ys.deriv,
                     "td-both" = alpha * Ys + Dalpha * Ys.deriv,
                     "shared-RE" = c(b %*% alpha))
        W <- object$x$W[object$id, seq_along(gammas), drop = FALSE]
        eta.tw <- if (ncol(W)) c(W %*% gammas) else rep(0, length(tt))
        cumHaz <- if (survMod == "weibull-PH") {
            Vi <- exp(log(sigma.t) + (sigma.t - 1) * log(c(t(st))) + tt)
            exp(eta.tw) * P * tapply(rep(wk, length.out = length(Vi)) * Vi, id.GK, sum)
        } else if (survMod == "spline-PH") {
            kn <- object$control$knots
            W2s <- splineDesign(unlist(kn, use.names = FALSE), c(t(st)), 
                                ord = object$control$ordSpline, outer.ok = TRUE)
            Vi <- exp(c(W2s %*% Bs.gammas) + tt)
            exp(eta.tw) * P * tapply(rep(wk, length.out = length(Vi)) * Vi, id.GK, sum)
        }
        names(cumHaz) <- names(object$y$y)
        cumHaz
    }
}
