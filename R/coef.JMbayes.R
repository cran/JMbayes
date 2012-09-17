coef.JMbayes <-
function (object, process = c("Longitudinal", "Event"), ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    process <- match.arg(process)
    if (process == "Longitudinal") {
        betas <- object$coefficients$betas
        out <- matrix(betas, nrow = length(object$y$logT), ncol = length(betas), byrow = TRUE)
        colnames(out) <- names(betas)
        EB <- object$ranef
        out[, colnames(EB)] <- out[, colnames(EB)] + EB
        out
    } else {
        gammas <- object$coefficients$gammas
        out <- c(gammas, "Assoct" = as.vector(object$coefficients$alphas), 
            "Assoct.s" = as.vector(object$coefficients$Dalphas))
        jj <- grep("Assoct[!^\\.s]", names(out))
        ii <- setdiff(grep("Assoct", names(out)), jj)
        if (length(ii) > 1) {
            nn <- names(object$coefficients$alpha)
            names(out)[ii] <- if (length(nn) == 1) "Assoct" else {
                if (nn[1] == "") 
                    c("Assoct", paste("Assoct", nn[-1], sep = ":"))
                else
                    paste("Assoct", nn, sep = ":")
            }
        }
        if (length(jj) > 1) {
            nn <- names(object$coefficients$Dalpha)
            names(out)[jj] <- if (length(nn) == 1) "Assoct.s" else {
                if (nn[1] == "") 
                    c("Assoct.s", paste("Assoct.s", nn[-1], sep = ":"))
                else
                    paste("Assoct.s", nn, sep = ":")
            }
        }
        if ((lag <- object$y$lag) > 0) {
            kk <- grep("Assoct", names(out), fixed = TRUE)
            names(out)[kk] <- paste(names(out)[kk], "(lag=", lag, ")", sep = "")
        }
        out
    }
}
