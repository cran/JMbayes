predict.JMbayes <-
function (object, newdata, type = c("Marginal", "Subject"),
    interval = c("none", "confidence", "prediction"), level = 0.95, idVar = "id", 
    FtTimes = NULL, M = 300, returnData = FALSE, scale = 1.6, 
    weight = rep(1, nrow(newdata)), ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    type <- match.arg(type)
    interval <- match.arg(interval)
    if (type == "Marginal") {
        TermsX <- delete.response(object$termsYx)
        mf <- model.frame(TermsX, data = newdata)
        form <- reformulate(attr(TermsX, "term.labels"))
        X <- model.matrix(form, data = mf)
        out <- c(X %*% object$coefficients$betas)
        names(out) <- row.names(newdata)
        if (interval == "prediction") {
            warning("\nfor type = 'Marginal' only confidence intervals are calculated.")
            interval <- "confidence"
        }
        if (interval == "confidence") {
            codaFit <- do.call(rbind, object$codaFit)
            ind <- grep("betas[", colnames(codaFit), fixed = TRUE)
            betas <- codaFit[, ind, drop = FALSE]
            preds <- X %*% t(betas)
            se.fit <- apply(preds, 1, sd)
            alpha <- 1 - level
            low <- apply(preds, 1, quantile, probs = alpha/2)
            up <- apply(preds, 1, quantile, probs = 1 - alpha/2)
            names(se.fit) <- names(low) <- names(up) <- row.names(newdata)
            out <- list(pred = out, se.fit = se.fit, low = low, upp = up)
        }
        if (returnData) {
            out <- if (is.list(out)) 
                cbind(newdata, do.call(cbind, out))
            else
                cbind(newdata, pred = out)
        }
    } else {
        if (is.null(newdata[[idVar]]))
            stop("'idVar' not in 'newdata.\n'")
        survMod <- object$survMod
        timeVar <- object$timeVar
        interFact <- object$interFact
        robust <- object$robust
        robust.b <- object$robust.b
        df <- object$df
        df.b <- object$df.b
        param <- object$param
        extraForm <- object$extraForm
        indFixed <- extraForm$indFixed
        indRandom <- extraForm$indRandom
        TermsX <- object$termsYx
        TermsZ <- object$termsYz
        TermsX.deriv <- object$termsYx.deriv
        TermsZ.deriv <- object$termsYz.deriv
        mfX <- model.frame(TermsX, data = newdata)
        mfZ <- model.frame(TermsZ, data = newdata)
        formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
        formYz <- object$formYz
        na.ind <- as.vector(attr(mfX, "na.action"))
        na.ind <- if (is.null(na.ind)) {
            rep(TRUE, nrow(newdata))
        } else {
            !seq_len(nrow(newdata)) %in% na.ind
        }
        id <- as.numeric(unclass(newdata[[idVar]]))
        id <- id. <- match(id, unique(id))
        id <- id[na.ind]
        y <- model.response(mfX)
        X <- model.matrix(formYx, mfX)
        Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
        TermsT <- object$termsT
        data.id <- newdata[!duplicated(id), ]
        idT <- data.id[[idVar]]
        idT <- match(idT, unique(idT))
        mfT <- model.frame(delete.response(TermsT), data = data.id)
        formT <- if (!is.null(kk <- attr(TermsT, "specials")$strata)) {
            strt <- eval(attr(TermsT, "variables"), data.id)[[kk]]
            tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
            reformulate(attr(tt, "term.labels"))
        } else if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
            tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
            reformulate(attr(tt, "term.labels"))
        } else {
            tt <- attr(delete.response(TermsT), "term.labels")
            if (length(tt)) reformulate(tt) else reformulate("1")
        }
        W <- model.matrix(formT, mfT)
        WintF.vl <- WintF.sl <- as.matrix(rep(1, nrow(data.id)))
        if (!is.null(interFact)) {
            if (!is.null(interFact$value))
                WintF.vl <- model.matrix(interFact$value, data = data.id)
            if (!is.null(interFact$slope))
                WintF.sl <- model.matrix(interFact$slope, data = data.id)
        }
        obs.times <- split(newdata[[timeVar]], id.)
        last.time <- tapply(newdata[[timeVar]], id., tail, n = 1)
        times.to.pred <- if (is.null(FtTimes)) {
            lapply(last.time, 
                function (t) seq(t, max(object$times) + 
                    0.1 * mad(object$times), length = 25))
        } else {
            if (!is.list(FtTimes) || length(FtTimes) != length(last.time))
                rep(list(FtTimes), length(last.time))
            else
                FtTimes
        }
        n <- length(object$y$logT)
        n.tp <- length(last.time)
        ncx <- ncol(X)
        ncz <- ncol(Z)
        ncww <- ncol(W)
        lag <- object$y$lag
        betas <- object$coefficients$betas
        sigma <- object$coefficients$sigma
        D <- object$coefficients$D
        diag.D <- ncol(D) == 1 & nrow(D) > 1
        D <- if (diag.D) diag(c(D)) else D
        gammas <- object$coefficients$gammas
        alpha <- object$coefficients$alpha
        Dalpha <- object$coefficients$Dalpha
        sigma.t <- object$coefficients$sigma.t
        Bs.gammas <- object$coefficients$Bs.gammas
        list.thetas <- list(betas = betas, sigma = sigma, 
            gammas = gammas, alpha = alpha, 
            Dalpha = Dalpha, sigma.t = sigma.t, 
            Bs.gammas = Bs.gammas, D = D)
        if (survMod  == "spline-PH") {
            if (ncww == 1) {
                W <- NULL
                ncww <- 0
            } else {
                W <- W[, -1, drop = FALSE]
                ncww <- ncww - 1
            }
        }
        list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
        thetas <- unlist(as.relistable(list.thetas))
        Var.thetas <- vcov(object)
        environment(log.posterior.b) <- environment(ModelMats) <- environment()
        # construct model matrices to calculate the survival functions
        obs.times.surv <- split(data.id[[timeVar]], idT)
        survMats.last <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            survMats.last[[i]] <- ModelMats(last.time[i], ii = i)
        }
        data.id2 <- newdata[!duplicated(id), ]
        data.id2 <- data.id2[rep(1:nrow(data.id2), 
            sapply(times.to.pred, length)), ]
        data.id2[[timeVar]] <- unlist(times.to.pred)
        mfXpred <- model.frame(TermsX, data = data.id2)
        mfZpred <- model.frame(TermsZ, data = data.id2)
        Xpred <- model.matrix(formYx, mfXpred)
        Zpred <- model.matrix(formYz, mfZpred)
        id2 <- as.numeric(unclass(data.id2[[idVar]]))
        id2 <- match(id2, unique(id2))
        # calculate the Empirical Bayes estimates and their (scaled) variance
        modes.b <- matrix(0, n.tp, ncz)
        Vars.b <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            betas.new <- betas
            sigma.new <- sigma
            D.new <- D
            gammas.new <- gammas
            alpha.new <- alpha
            Dalpha.new <- Dalpha
            sigma.t.new <- sigma.t
            Bs.gammas.new <- Bs.gammas
            ff <- function (b, y, tt, mm, i) 
                -log.posterior.b(b, y, Mats = tt, survMod = mm, ii = i)
            opt <- try(optim(rep(0, ncz), ff, y = y, tt = survMats.last, mm = survMod, i = i, 
                method = "BFGS", hessian = TRUE), TRUE)
            if (inherits(opt, "try-error")) {
                gg <- function (b, y, tt, mm, i) cd(b, ff, y = y, tt = tt, mm = mm, i = i)
                opt <- optim(rep(0, ncz), ff, gg, y = y, tt = survMats.last, mm = survMod, 
                    i = i, method = "BFGS", hessian = TRUE)
            }
            modes.b[i, ] <- opt$par
            Vars.b[[i]] <- scale * solve(opt$hessian)        
        }
        res <- vector("list", M)
        success.rate <- matrix(FALSE, M, n.tp)
        b.old <- b.new <- modes.b
        if (n.tp == 1)
            dim(b.old) <- dim(b.new) <- c(1, ncz)    
        codaFit <- do.call("rbind", object$codaFit)
        if (M > nrow(codaFit)) {
            warning("\n'M' cannot be set as greater than ", nrow(codaFit))
            M <- nrow(codaFit)
            res <- vector("list", M)
            success.rate <- matrix(FALSE, M, n.tp)
        }
        codaFit <- codaFit[sample(nrow(codaFit), M), , drop = FALSE]
        codaFit <- cbind(codaFit, "sigma" = 1 / codaFit[, "tau"])
        ind.thetas <- sapply(names(list.thetas), grep, x = colnames(codaFit), fixed = TRUE)
        if (!is.null(ind.thetas$Bs.gammas) && !is.null(object$coefficients$gammas))
            ind.thetas$gammas <- setdiff(ind.thetas$gammas, ind.thetas$Bs.gammas)
        if (param %in% c("td-extra", "td-both"))
            ind.thetas$D <- setdiff(ind.thetas$D, ind.thetas$Dalpha)
        if (param == "td-both")
            ind.thetas$alpha <- setdiff(ind.thetas$alpha, ind.thetas$Dalpha)
        if (param == "td-extra")
            ind.thetas$alpha <- NULL
        if (param %in% c("td-value", "shared-RE"))
            ind.thetas$Dalpha <- NULL
        if (survMod == "weibull-PH") {
            ind.thetas$sigma <- setdiff(ind.thetas$sigma, ind.thetas$sigma.t)
            ind.thetas$Bs.gammas <- NULL
        }
        for (m in 1:M) {
            # Step 1: simulate new parameter values
            thetas.new <- lapply(ind.thetas, function (k) codaFit[m, k])
            thetas.new$D <- solve(matrix(thetas.new$D, ncz, ncz, TRUE))
            betas.new <- thetas.new$betas
            sigma.new <- thetas.new$sigma
            gammas.new <- thetas.new$gammas
            alpha.new <- thetas.new$alpha
            Dalpha.new <- thetas.new$Dalpha
            D.new <- thetas.new$D
            sigma.t.new <- thetas.new$sigma.t
            Bs.gammas.new <- thetas.new$Bs.gammas
            y.new <- vector("list", n.tp)
            for (i in seq_len(n.tp)) {
                # Step 2: simulate new random effects values
                proposed.b <- rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
                dmvt.old <- dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
                dmvt.proposed <- dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
                a <- min(exp(log.posterior.b(proposed.b, y, survMats.last, survMod, ii = i) + dmvt.old - 
                        log.posterior.b(b.old[i, ], y, survMats.last, survMod, ii = i) - dmvt.proposed), 1)
                ind <- runif(1) <= a
                success.rate[m, i] <- ind
                if (!is.na(ind) && ind)
                    b.new[i, ] <- proposed.b
                # Step 3: compute future Ys
                Xpred.i <- Xpred[id2 == i, , drop = FALSE]
                Zpred.i <- Zpred[id2 == i, , drop = FALSE]
                mu.i <- as.vector(c(Xpred.i %*% betas.new) + 
                    rowSums(Zpred.i * rep(b.new[i, ], each = nrow(Zpred.i))))
                y.new[[i]] <- if (interval == "confidence") weight[i] * mu.i else 
                    if (interval == "prediction") {
                        if (!robust) 
                            weight[i] * rnorm(length(mu.i), mu.i, sigma.new)
                        else
                            weight[i] * rgt(length(mu.i), mu.i, sigma.new, df)
                    }
            }
            b.old <- b.new
            res[[m]] <- y.new
        }
        oo <- vector("list", n.tp)
        for (i in seq_len(n.tp)) {
            oo[[i]] <- do.call(rbind, sapply(res, "[", i))
        }
        out <- weight[i] * as.vector(c(Xpred %*% betas) + 
            rowSums(Zpred * modes.b[id2, , drop = FALSE]))
        if (interval %in% c("confidence", "prediction")) {
            alpha <- 1 - level
            se.fit <- lapply(oo, function (m) {
                if (is.matrix(m)) 
                    apply(m, 2, sd)
                else 
                    sd(m)
            })
            f1 <- function (mat) apply(mat, 2, quantile, probs = alpha/2)
            f2 <- function (mat) apply(mat, 2, quantile, probs = 1 - alpha/2)
            low <- lapply(oo, f1) 
            up <- lapply(oo, f2)
            out <- list(pred = out, se.fit = unlist(se.fit), 
                low = unlist(low), upp = unlist(up))
        }
        if (returnData) {
            newdata$pred <- c(X %*% betas) + rowSums(Z * modes.b[id, ])
            out <- if (is.list(out)) {
                newdata$upp <- newdata$low <- newdata$se.fit <- NA
                rbind(newdata, cbind(data.id2, do.call(cbind, out)))
            } else {
                rbind(newdata, cbind(data.id2, pred = out))
            }
        } else
            attr(out, "time.to.pred") <- times.to.pred
    }
    class(out) <- c(class(out), "predict.JMbayes")
    out
}
