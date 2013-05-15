print.summary.JMbayes <-
function (x, digits = max(4, getOption("digits") - 4), 
        printKnots = FALSE, ...) {
    if (!inherits(x, "summary.JMbayes"))
        stop("Use only with 'summary.JMbayes' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", 
        collapse = "\n"), "\n\n", sep = "")
    cat("Data Descriptives:\n")
    pcEv <- round(100 * sum(x$d) / x$n, 1)
    cat("Longitudinal Process\t\tEvent Process")
    cat("\nNumber of Observations: ", x$N, "\tNumber of Events: ", 
        sum(x$d), " (", pcEv, "%)", sep = "")
    cat("\nNumber of Groups:", length(unique(x$id)))
    cat("\n\nJoint Model Summary:")
    if (!x$robust && !x$robust.b) {
        cat("\nLongitudinal Process: Linear mixed-effects model")
    } else if (x$robust && !x$robust.b) {
        cat("\nLongitudinal Process: Linear mixed-effects model with Student's-t(df=", x$df, ") errors", sep = "")
    } else if (!x$robust && x$robust.b) {
        cat("\nLongitudinal Process: Linear mixed-effects model with Student's-t(df=", x$df.b, ") random effects", sep = "")
    } else {
        cat("\nLongitudinal Process: Linear mixed-effects model with Student's-t(df=", x$df, ") errors and\n\t\t",
            "Student's-t(df=", x$df.b, ") random effects", sep = "")
    }
    cat("\nEvent Process: ")
    if (x$survMod == "weibull-PH") {
        cat("Weibull relative risk model\n")
    } else if (x$survMod == "spline-PH") {
        xx <- if (length(x$control$knots) == 1) {
            kk <- round(unique(x$control$knots[[1]]), 1)
            paste(kk[-c(1, length(kk))], collapse = ", ")
        } else {
            paste(names(x$control$knots), sapply(x$control$knots, function (k) {
                kk <- round(unique(k), 1)
                paste(kk[-c(1, length(kk))], collapse = ", ")
            }), sep = ": ", collapse = "\n\t\t")
        }
        if (printKnots)
            cat("Relative risk model with spline-approximated baseline risk function (knots at: ", xx, ")\n", sep = "")
        else
            cat("Relative risk model with spline-approximated\n\t\tbaseline risk function\n")
    }
    cat("Parameterization:", switch(x$param, "td-value" = "Time-dependent", 
        "td-extra" = "Time-dependent slope", "td-both" = "Time-dependent + time-dependent slope",
        "shared-RE" = "shared random effects"), "\n\n")
    if (!is.null(x$DIC)){ 
        model.sum <- data.frame(logLik = x$logLik, DIC = x$DIC, pD = x$pD, row.names = "")
        print(model.sum)
    }
    cat("\nVariance Components:\n")
    D <- x$D
    ncz <- nrow(D)
    diag.D <- ncz != ncol(D)
    sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- as.data.frame(round(rbind(sds, "Residual" = x$sigma), digits))
            names(dat) <- "StdDev"
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- rbind(mat, c(x$sigma, rep(0, ncz - 1)))
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(D)[-ncz], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- c(dimnames(D)[[1]], "Residual")
        }
    } else {
        dat <- data.frame("StdDev" = c(sds, x$sigma), row.names = c(rownames(D), "Residual"), 
            check.rows = FALSE, check.names = FALSE)
    }
    print(dat)
    cat("\nCoefficients:")
    cat("\nLongitudinal Process\n")
    out <- as.data.frame(round(x$"CoefTable-Long", digits))
    print(out)
    cat("\nEvent Process\n")
    out <- as.data.frame(round(x$"CoefTable-Event", digits))
    print(out)
    cat("\nMCMC summary:\n")
    cat("program:", switch(x$control$program, "WinBUGS" =, "winbugs" = "WinBUGS",
        "OpenBUGS" =, "openbugs" = "OpenBUGS", "JAGS" =, "jags" = "JAGS"), "\n")
    cat("chains:", x$control$n.chains, "\nthinning:", x$control$n.thin, "\n")
    cat("iterations:", x$control$n.iter, "\nburn-in:", x$control$n.burnin)
    cat("\n")
    invisible(x)    
}
