summary.JMbayes <-
function (object, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    coefs <- object$coefficients
    sds <- object$StErr
    CIs <- object$CIs
    coefsY <- cbind("Value" = coefs$betas, "Std.Err" = sds$betas, "2.5%" = CIs$betas[1, ], 
        "97.5%" = CIs$betas[2, ])
    gammas <- if (object$survMod == "weibull-PH") {
        c(coefs$gammas, "Assoct" = as.vector(coefs$alpha),
            "Assoct.s" = as.vector(coefs$Dalpha), 
            "shape" = as.vector(coefs$sigma.t))
    } else {
        c(coefs$gammas, "Assoct" = as.vector(coefs$alpha),
            "Assoct.s" = as.vector(coefs$Dalpha), 
            "Bs.gammas" = as.vector(coefs$Bs.gammas))
    }
    if (object$param == "shared-RE") {
        ii <- grep("Assoct", names(gammas), fixed = TRUE)
        names(gammas)[ii] <- paste("Assoct.", colnames(object$x$Z), sep = "")
    }
    sds.gammas <- c(sds$gammas, sds$alpha, sds$Dalpha, sds$sigma.t, sds$Bs.gammas)
    cis.gammas <- rbind(c(CIs$gammas[1, ], if (is.matrix(CIs$alpha)) CIs$alpha[1, ] else CIs$alpha[1], 
            CIs$Dalpha[1], CIs$sigma.t[1], CIs$Bs.gammas[1, ]), 
        c(CIs$gammas[2, ], if (is.matrix(CIs$alpha)) CIs$alpha[2, ] else CIs$alpha[2], 
            CIs$Dalpha[2], CIs$sigma.t[2], CIs$Bs.gammas[2, ]))
    coefsT <- cbind("Value" = gammas, "Std.Err" = sds.gammas, "2.5%" = cis.gammas[1, ],
        "97.5%" = cis.gammas[2, ])
    out <- list("CoefTable-Long" = coefsY, "CoefTable-Event" = coefsT, 
        D = coefs$D, sigma = coefs$sigma)
    out$N <- nrow(object$x$X)
    out$n <- length(object$y$logT)
    out$d <- object$y$event
    out$id <- object$id
    out$survMod <- object$survMod
    out$control <- object$control
    out$knots <- unique(object$knots)
    out$conv <- object$conv
    out$param <- object$param
    out$robust <- object$robust
    out$df <- object$df
    out$DIC <- object$DIC
    out$pD <- object$pD
    out$logLik <- logLik(object)
    out$call <- object$call
    class(out) <- "summary.JMbayes"
    out
}
