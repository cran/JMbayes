jointModelBayes <-
function (lmeObject, survObject, timeVar, survMod = c("weibull-PH", "spline-PH"), 
        param = c("td-value", "td-extra", "td-both", "shared-RE"), robust = FALSE, robust.b = FALSE, 
        df = 4, df.b = 4, lag = 0, init = NULL, extraForm = NULL, priors = NULL, control = list(), ...) {
    cl <- match.call()
    if (!inherits(lmeObject, "lme"))
        stop("\n'lmeObject' must inherit from class lme.")
    if (length(lmeObject$group) > 1)
        stop("\nnested random-effects are not allowed in lme().")
    if (!is.null(lmeObject$modelStruct$corStruct))
        warning("correlation structure in 'lmeObject' is ignored.\n")
    if (!is.null(lmeObject$modelStruct$varStruct))
        warning("variance structure in 'lmeObject' is ignored.\n")        
    if (!inherits(survObject, "coxph") && !inherits(survObject, "survreg"))
        stop("\n'survObject' must inherit from class coxph or class survreg.")
    if (is.null(survObject$x))
        stop("\nuse argument 'x = TRUE' in ", if (inherits(survObject, "coxph")) "'coxph()'." else "'survreg()'.")
    if (length(timeVar) != 1 || !is.character(timeVar))
        stop("\n'timeVar' must be a character string.")
    if (param %in% c("td-extra", "td-both") && is.null(extraForm)) {
        stop("\nwhen parameterization is 'td-extra' or 'td-both' you need to specify the 'extraForm' argument.")
    }
    if (!param %in% c("td-extra", "td-both") && !is.null(extraForm)) {
        stop("\nyou have defined 'extraForm' but the parameterization is neither 'td-extra' nor 'td-both'.")
    }    
    if (param %in% c("td-extra", "td-both") && !is.list(extraForm)) {
        stop("\nthe 'extraForm' argument must be a list with components 'fixed' (a formula),\n\t'indFixed'", 
            "(a numeric vector), 'random' (a formula) and 'indRandom' (a numeric vector).")
    }
    survMod <- match.arg(survMod)
    param <- match.arg(param)
    # extract response & design matrix survival process
    formT <- formula(survObject)
    if (inherits(survObject, "coxph")) {
        W <- survObject$x
        Time <- survObject$y[, 1]
    } else {
        W <- survObject$x[, -1, drop = FALSE]
        Time <- exp(survObject$y[, 1])
    }
    nT <- length(Time)
    if (!length(W))
        W <- NULL
    event <- survObject$y[, 2]
    zeros <- numeric(nT) # for the zeros trick
    # longitudinal process
    id <- as.vector(unclass(lmeObject$groups[[1]]))
    b <- data.matrix(ranef(lmeObject))
    one.RE <- ncol(b) == 1
    if (one.RE)
        b <- cbind(b, rep(0, nrow(b)))
    dimnames(b) <- NULL
    nY <- nrow(b)
    if (nY != nT)
        stop("sample sizes in the longitudinal and event processes differ.\n")
    TermsX <- lmeObject$terms
    data <- lmeObject$data[all.vars(TermsX)]
    data <- data[complete.cases(data), ]
    formYx <- formula(lmeObject)
    mfX <- model.frame(TermsX, data = data)
    X <- model.matrix(formYx, mfX)
    formYz <- formula(lmeObject$modelStruct$reStruct[[1]])    
    mfZ <- model.frame(terms(formYz), data = data)
    TermsZ <- attr(mfZ, "terms")
    Z <- model.matrix(formYz, mfZ)
    if (one.RE)
        Z <- cbind(Z, rep(0, nrow(Z)))
    y.long <- model.response(mfX, "numeric")
    offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
    # control values
    con <- list(program = "JAGS", n.chains = 1, n.iter = 10000, n.burnin = 5000, 
        n.thin = 1, n.adapt = 1000, K = 100, C = 5000, working.directory = getwd(), 
        bugs.directory = "C:/Program Files/WinBUGS14/", openbugs.directory = NULL, 
        clearWD = TRUE, over.relax = TRUE, knots = NULL, ObsTimes.knots = TRUE, 
        lng.in.kn = 5, ordSpline = 4, bugs.seed = 1, quiet = FALSE)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    data.id <- data[!duplicated(id), ]
    if (!timeVar %in% names(data))
        stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'lmeObject'.")
    # extra design matrices for the longitudinal part
    data.id[[timeVar]] <- pmax(Time - lag, 0)
    if (param %in% c("td-value", "td-both")) {
        mfX.id <- model.frame(TermsX, data = data.id)
        mfZ.id <- model.frame(TermsZ, data = data.id)
        Xtime <- model.matrix(formYx, mfX.id)
        Ztime <- model.matrix(formYz, mfZ.id)
        if (one.RE)
            Ztime <- cbind(Ztime, rep(0, nrow(Ztime)))
        long <- c(X %*% fixef(lmeObject)) + rowSums(Z * b[id, ])
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.deriv <- model.frame(terms(extraForm$fixed), data = data)
        TermsX.deriv <- attr(mfX.deriv, "terms")
        mfZ.deriv <- model.frame(terms(extraForm$random), data = data)
        TermsZ.deriv <- attr(mfZ.deriv, "terms")
        mfX.deriv.id <- model.frame(TermsX.deriv, data = data.id)
        mfZ.deriv.id <- model.frame(TermsZ.deriv, data = data.id)      
        Xtime.deriv <- model.matrix(extraForm$fixed, mfX.deriv.id)
        Ztime.deriv <- model.matrix(extraForm$random, mfZ.deriv.id)
        Xderiv <- model.matrix(extraForm$fixed, mfX.deriv)
        Zderiv <- model.matrix(extraForm$random, mfZ.deriv)        
        long.deriv <- as.vector(c(Xderiv %*% fixef(lmeObject)[extraForm$indFixed]) + 
            if (length(extraForm$indRandom) > 1 || extraForm$indRandom) 
                rowSums(Zderiv * b[id, extraForm$indRandom, drop = FALSE])
            else
                rep(0, nrow(Zderiv)))
    }
    if (param == "td-value")
        long.deriv <- NULL
    if (param == "td-extra")
        long <- NULL
    # response vectors and design matrices
    y <- list(y = y.long, offset = offset, logT = log(Time), event = event, zeros = zeros, lag = lag)
    x <- list(X = X, Z = Z, 
        W = if (survMod == "weibull-PH") { 
            if (is.null(W)) cbind(rep(1, nT), rep(0, nT)) else cbind(1, W)
        } else {
            if (is.null(W)) cbind(rep(0, nT), rep(0, nT)) else {
                if (ncol(W) == 1) cbind(W, rep(0, nT)) else W
            }
    })
    x <- switch(param, 
        "td-value" = c(x, list(Xtime = Xtime, Ztime = Ztime)),
        "td-extra" = c(x, list(Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)),
        "td-both" = c(x, list(Xtime = Xtime, Ztime = Ztime, 
            Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)),
        "shared-RE" = x)
    wk <- gaussKronrod()$wk
    sk <- gaussKronrod()$sk
    K <- length(sk)
    P <- Time/2
    st <- outer(P, sk + 1)
    id.GK <- rep(seq_along(Time), each = K)
    data.id2 <- data.id[id.GK, ]
    data.id2[[timeVar]] <- c(t(st))
    x <- c(x, list(P = P, st = st, wk = wk))    
    if (param %in% c("td-value", "td-both")) {
        mfX <- model.frame(TermsX, data = data.id2)
        mfZ <- model.frame(TermsZ, data = data.id2)
        Xs <- model.matrix(formYx, mfX)
        Zs <- model.matrix(formYz, mfZ)
        if (one.RE)
            Zs <- cbind(Zs, rep(0, nrow(Zs)))
        x <- c(x, list(Xs = Xs, Zs = Zs))
    }
    if (param %in% c("td-extra", "td-both")) {
        mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
        mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
        Xs.deriv <- model.matrix(extraForm$fixed, mfX.deriv)
        Zs.deriv <- model.matrix(extraForm$random, mfZ.deriv)
        x <- c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv))
    }
    evTime <- cens <- Time
    cens[event == 1] <- 0
    evTime[event == 0] <- NA
    y <- c(y, list(evTime = evTime, censTime = cens))    
    # extra design matrices for the log approximated baseline hazard
    if (survMod == "spline-PH") { 
        kn <- if (is.null(con$knots)) {
            pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
            pp <- tail(head(pp, -1), -1)
            tt <- if (con$ObsTimes.knots) Time else Time[event == 1]
            quantile(tt, pp, names = FALSE)
        } else {
            con$knots
        }
        kn <- kn[kn < max(Time)]
        rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
        con$knots <- rr
        W2 <- splineDesign(rr, Time, ord = con$ordSpline)
        if (any(colSums(W2) == 0))
            stop("\nsome of the knots of the B-splines basis are set outside the range",
                "\n   of the observed event times for one of the strata; refit the model", 
                "\n   setting the control argument 'equal.strata.knots' to FALSE.")
        W2s <- splineDesign(rr, c(t(st)), ord = con$ordSpline)
        x <- c(x, list(W2 = W2, W2s = W2s))
    }
    # default priors
    ncX <- ncol(X)
    ncZ <- ncol(Z)
    ncW <- ncol(x$W)
    ncW2 <- ncol(x$W2)
    betas <- rep(0, ncX)
    var.betas <- rep(con$K, ncX)
    sigma2 <- lmeObject$sigma^2
    VC <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*", 
        lmeObject$sigma^2)[[1]]
    if (one.RE)
        VC <- cbind(c(VC, 0), c(0, 1))
    inv.VC <- solve(VC)
    alphas <- Dalphas <- 0
    var.alphas <- var.Dalphas <- con$K
    gammas <- rep(0, ncW)
    var.gammas <- rep(con$K, ncW)
    if (survMod == "weibull-PH") {
        if (is.null(W)) var.gammas[2] <- 1e-05
    } else {
        if (is.null(W)) var.gammas[1:2] <- 1e-05
        if (!is.null(W) && ncol(W) == 1) var.gammas[2] <- 1e-05
        Bs.gammas <- rep(0, ncW2)
        var.Bs.gammas <- rep(con$K/10, ncW2)
    }
    # Data to be passed to WinBUGS
    Data <- list(N = nY, df = df, df.b = df.b, C = con$C, K = K, R = con$R, ncX = ncX, 
        ncZ = ncZ, ncW = ncW, ncW2 = ncW2,
        priorMean.betas = betas, priorTau.betas = diag(1 / var.betas),
        priorA.tau = (1/sigma2)^2 / 10, priorB.tau = (1/sigma2) / 10,
        priorMean.gammas = gammas, priorTau.gammas = diag(1 / var.gammas),
        priorMean.alphas = alphas, priorTau.alphas = 1 / var.alphas,
        priorMean.Dalphas = alphas, priorTau.Dalphas = 1 / var.Dalphas,
        priorA.sigma.t = 10, priorB.sigma.t = 10, priorA.rho = 2, priorB.rho = 1,
        nb = ncZ, mu0 = rep(0, ncZ), priorR.D = ncZ * inv.VC, priorK.D = ncZ)
    Data <- c(Data, x, y)
    if (param %in% c("td-extra", "td-both") && !is.null(extraForm))
        Data <- c(Data, list(indFixed = extraForm$indFixed, 
            indRandom = extraForm$indRandom,
            ncX.deriv = ncol(Xtime.deriv), ncZ.deriv = ncol(Ztime.deriv)))
    if (survMod == "weibull-PH") {
        Data <- c(Data, list(Time = Time, logTime = log(Time)))
    } else {
        Data <- c(Data, list(priorMean.Bs.gammas = Bs.gammas, 
            priorTau.Bs.gammas = diag(1 / var.Bs.gammas)))
    }
    if (param == "shared-RE") {
        Data$priorMean.alphas <- rep(0, ncZ)
        Data$priorTau.alphas <- diag(1 / rep(con$K, ncZ))
    }
    # user priors
    if (!is.null(priors)) {
        if (!is.list(priors) || !names(priors) %in% names(Data)) {
            warning("'priors' is not a list with appropriate names (check the help file); ", 
                "default priors are used instead.\n")
        } else { 
            Data[names(priors)] <- priors
        }
    }
    # initial values
    initial.values <- list(betas = fixef(lmeObject), tau = 1/sigma2, 
        inv.D = inv.VC, gammas = gammas, b = b)
    initial.values$alphas <- initial.values$Dalphas <- 0
    if (param == "shared-RE")
        initial.values$alphas <- rep(0, ncZ)
    if (survMod == "weibull-PH") {
        initial.values$sigma.t <- 1
    } else {
        initial.values$Bs.gammas <- rep(0.1, ncW2)
    }
    if (!is.null(init)) {
        lngths <- lapply(initial.values[(nam.init <- names(init))], length)
        if (!is.list(init) || !isTRUE(all.equal(lngths, lapply(init, length)))) {
            warning("'init' is not a list with elements numeric vectors of appropriate ", 
                "length; default starting values are used instead.\n")
        } else {
            initial.values[nam.init] <- init
        }
    }
    # joint model fit
    model <- switch(paste(survMod, param, robust, robust.b, sep = "/"),
        "weibull-PH/td-value/FALSE/FALSE" = Weib.td.value,
        "weibull-PH/td-extra/FALSE/FALSE" = Weib.td.slope,
        "weibull-PH/td-both/FALSE/FALSE" = Weib.td.both,
        "weibull-PH/shared-RE/FALSE/FALSE" = Weib.sharedRE,
        "spline-PH/td-value/FALSE/FALSE" = spline.td.value,
        "spline-PH/td-extra/FALSE/FALSE" = spline.td.slope,
        "spline-PH/td-both/FALSE/FALSE" = spline.td.both,
        "spline-PH/shared-RE/FALSE/FALSE" = spline.sharedRE,
        ###            
        "weibull-PH/td-value/TRUE/FALSE" = WeibRob.y.td.value,
        "weibull-PH/td-extra/TRUE/FALSE" = WeibRob.y.td.slope,
        "weibull-PH/td-both/TRUE/FALSE" = WeibRob.y.td.both,
        "weibull-PH/shared-RE/TRUE/FALSE" = WeibRob.y.sharedRE,
        "spline-PH/td-value/TRUE/FALSE" = splineRob.y.td.value,
        "spline-PH/td-extra/TRUE/FALSE" = splineRob.y.td.slope,
        "spline-PH/td-both/TRUE/FALSE" = splineRob.y.td.both,
        "spline-PH/shared-RE/TRUE/FALSE" = splineRob.y.sharedRE,
        ###
        "weibull-PH/td-value/FALSE/TRUE" = WeibRob.b.td.value,
        "weibull-PH/td-extra/FALSE/TRUE" = WeibRob.b.td.slope,
        "weibull-PH/td-both/FALSE/TRUE" = WeibRob.b.td.both,
        "weibull-PH/shared-RE/FALSE/TRUE" = WeibRob.b.sharedRE,
        "spline-PH/td-value/FALSE/TRUE" = splineRob.b.td.value,
        "spline-PH/td-extra/FALSE/TRUE" = splineRob.b.td.slope,
        "spline-PH/td-both/FALSE/TRUE" = splineRob.b.td.both,
        "spline-PH/shared-RE/FALSE/TRUE" = splineRob.b.sharedRE,
        ###            
        "weibull-PH/td-value/TRUE/TRUE" = WeibRob.yb.td.value,
        "weibull-PH/td-extra/TRUE/TRUE" = WeibRob.yb.td.slope,
        "weibull-PH/td-both/TRUE/TRUE" = WeibRob.yb.td.both,
        "weibull-PH/shared-RE/TRUE/TRUE" = WeibRob.yb.sharedRE,
        "spline-PH/td-value/TRUE/TRUE" = splineRob.yb.td.value,
        "spline-PH/td-extra/TRUE/TRUE" = splineRob.yb.td.slope,
        "spline-PH/td-both/TRUE/TRUE" = splineRob.yb.td.both,
        "spline-PH/shared-RE/TRUE/TRUE" = splineRob.yb.sharedRE)
    namesModel <- namesJMbayes(survMod, param, robust, robust.b)
    parms <- c("betas", "tau", "inv.D", "gammas", "b")
    parms <- switch(paste(survMod, param, sep = "/"),
        "weibull-PH/td-value" = c(parms, "alphas", "sigma.t"),
        "weibull-PH/td-extra" = c(parms, "Dalphas", "sigma.t"),
        "weibull-PH/td-both" = c(parms, "alphas", "Dalphas", "sigma.t"),
        "weibull-PH/shared-RE" = c(parms, "alphas", "sigma.t"),
        "spline-PH/td-value" = c(parms, "alphas", "Bs.gammas"),
        "spline-PH/td-extra" = c(parms, "Dalphas", "Bs.gammas"),
        "spline-PH/td-both" = c(parms, "alphas", "Dalphas", "Bs.gammas"),
        "spline-PH/shared-RE" = c(parms, "alphas", "Bs.gammas"))
    write.model.JMbayes(model, file.path(con$working.directory, "JM.txt"), Data, 
        param, extraForm, con$program)
    Data <- Data[names(Data) %in% namesModel]
    initial.values <- initial.values[names(initial.values) %in% parms]
    fit <- if (con$program %in% c("WinBUGS", "winbugs")) {
        if (!require("R2WinBUGS")) 
            stop("'R2WinBUGS' is required.\n")
        R2WinBUGS::bugs(Data, list(initial.values), parms, model.file = "JM.txt", 
            n.chains = con$n.chains, n.iter = con$n.iter, n.burnin = con$n.burnin, 
            n.thin = con$n.thin, working.directory = con$working.directory, 
            clearWD = con$clearWD, bugs.directory = con$bugs.directory, 
            over.relax = con$over.relax, bugs.seed = con$bugs.seed)
    } else if (con$program %in% c("OpenBUGS", "openbugs")) {
        if (!require("R2OpenBUGS")) 
            stop("'R2OpenBUGS' is required.\n")
        R2OpenBUGS::bugs(Data, list(initial.values), parms, model.file = "JM.txt", 
            n.chains = con$n.chains, n.iter = con$n.iter, n.burnin = con$n.burnin, 
            n.thin = con$n.thin, clearWD = con$clearWD, over.relax = con$over.relax, 
            bugs.seed = con$bugs.seed, working.directory = con$working.directory,
            OpenBUGS.pgm = con$openbugs.directory)
    } else if (con$program %in% c("JAGS", "jags")) {
        if (!require("rjags")) 
            stop("'rjags' is required.\n")
        JMjags.model <- jags.model(file = "JM.txt", data = Data, inits = list(initial.values),
            n.chains = con$n.chains, n.adapt = con$n.adapt, quiet = con$quiet)
        update(JMjags.model, con$n.burnin)
        coda.samples(JMjags.model, parms, n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
    } else {
        stop("'program' should be one of 'WinBUGS', 'winbugs', 'OpenBUGS', 'openbugs', ",
            "'JAGS', 'jags'.\n")
    }
    file.remove(file.path(con$working.directory, "JM.txt"))
    # output result
    codaFit <- as.mcmc.list(fit)
    Bs <- do.call(rbind, codaFit)
    ind.bs <- grep("b[", colnames(Bs), fixed = TRUE)
    ranef <- Bs[, ind.bs, drop = FALSE]
    ord.col <- sapply(strsplit(colnames(ranef), "[", fixed = TRUE), "[", 2)
    ord.col <- sapply(strsplit(ord.col, ",", fixed = TRUE), "[", 1)
    ord.col <- order(as.numeric(ord.col))
    ranef <- ranef[, ord.col, drop = FALSE]
    postMeans <- matrix(colMeans(ranef), ncol = ncol(Z), byrow = TRUE)
    dimnames(postMeans) <- list(levels(factor(lmeObject$groups[[1]])), colnames(Z))
    postVars <- vector("list", nrow(postMeans))
    ind.var <- matrix(seq_len(ncol(ranef)), ncol = ncol(Z), byrow = TRUE)
    for (i in seq_along(postVars))
        postVars[[i]] <- var(ranef[, ind.var[i, ]])
    postVars[] <- lapply(postVars, function (m) { 
        dimnames(m) <- list(colnames(postMeans), colnames(postMeans))
        m
    })
    names(postVars) <- rownames(postMeans)
    if (one.RE) {
        postMeans <- postMeans[, 1, drop = FALSE]
        postVars <- lapply(postVars, function (m) m[1, 1, drop = FALSE])
    }
    #
    n.sims <- nrow(Bs)
    sims.list <- vector("list", length(parms))
    names(sims.list) <- parms
    for (p in seq_along(parms)) {
        ii <- grep(paste("^", parms[p], sep = ""), colnames(Bs))
        sims.list[[p]] <- Bs[, ii]
    }
    sims.list$b <- NULL
    h <- function (b, thetas) {
        thetas <- relist(thetas, skeleton = list.thetas)
        betas <- thetas$betas
        sigma <- thetas$sigma
        gammas <- thetas$gammas
        alphas <- thetas$alphas
        Dalphas <- thetas$Dalphas
        sigma.t <- thetas$sigma.t
        D <- thetas$D
        n <- length(y$logT)
        mu.y <- c(x$X %*% betas) + rowSums(x$Z * b[id, , drop = FALSE])
        log.p.y.b <- if (!robust) dnorm(y$y, mu.y, sigma, log = TRUE) 
            else dgt(y$y, mu.y, sigma, df, log = TRUE)
        log.p.y.b <- tapply(log.p.y.b, id, sum)
        log.p.b <- if (!robust.b) dmvnorm(b, rep(0, ncol(b)), D, log = TRUE)
            else dmvt(b, rep(0, ncol(b)), D, df.b, log = TRUE)
        eta.t <- if (!is.null(gammas)) c(x$W %*% rep(gammas, length.out = ncol(x$W))) 
            else rep(0, n)
        id.GK <- rep(seq_along(y$logT), each = 15)
        wk.long <- rep(wk, n)
        if (param %in% c("td-value", "td-both")) {
            Y <- c(x$Xtime %*% betas) + rowSums(x$Ztime * b)
            Ys <- c(x$Xs %*% betas) + rowSums(x$Zs * b[id.GK, , drop = FALSE])
        }
        if (param %in% c("td-extra", "td-both")) {
            Yderiv <- c(x$Xtime.deriv %*% betas[extraForm$indFixed]) + 
                rowSums(x$Ztime.deriv * b[, extraForm$indRandom, drop = FALSE])
            Ys.deriv <- c(x$Xs.deriv %*% betas[extraForm$indFixed]) + 
                rowSums(x$Zs.deriv * b[id.GK, extraForm$indRandom, drop = FALSE])
        }
        longSurv <- switch(param,
            "td-value" = alphas * Y, 
            "td-extra" = Dalphas * Yderiv,
            "td-both" = alphas * Y + Dalphas * Yderiv,
            "shared-RE" = c(b %*% alphas))
        longSurv.s <- switch(param,
            "td-value" = alphas * Ys, 
            "td-extra" = Dalphas * Ys.deriv,
            "td-both" = alphas * Ys + Dalphas * Ys.deriv,
            "shared-RE" = c(b %*% alphas)[id.GK])
        if (survMod == "weibull-PH") {
            log.hazard <- log(sigma.t) + (sigma.t - 1) * y$logT + eta.t + longSurv
            log.survival <- - exp(eta.t) * P * tapply(wk.long * exp(log(sigma.t) + 
                (sigma.t - 1) * log(c(t(st))) + longSurv.s), id.GK, sum)
        } else {
            log.hazard <- c(x$W2 %*% Bs.gammas) + eta.t + longSurv
            log.survival <- - exp(eta.t) * P * tapply(wk.long * exp(c(W2s %*% Bs.gammas) + 
                longSurv.s), id.GK, sum)
        }
        log.p.t.b <- event * log.hazard + log.survival
        - 2 * sum(log.p.y.b + log.p.t.b + log.p.b, na.rm = TRUE)
    }
    out <- list(codaFit = lapply(codaFit, function (x) x[, -ind.bs, drop = FALSE]),
        postMeans = postMeans, postVars = postVars)
    class(out$codaFit) <- "mcmc.list"
    ncz <- ncol(Z)
    indD <- cbind(rep(1:ncz, each = ncz), rep(1:ncz, ncz))
    DD <- t(sapply(seq_len(n.sims), function (i) {
        m <- matrix(0, ncz, ncz)
        m[indD] <- sims.list$inv.D[i, ]
        d <- solve(m)
        d[lower.tri(d, TRUE)]
    }))
    sims.list <- list(betas = sims.list$betas, sigma = sqrt(1 / sims.list$tau),
        gammas = sims.list$gammas, alphas = sims.list$alphas, 
        Dalphas = sims.list$Dalphas, sigma.t = sims.list$sigma.t,
        Bs.gammas = sims.list$Bs.gammas, D = DD)
    if (is.null(W)) {
        sims.list$gammas <- if (survMod == "weibull-PH") sims.list$gammas[, 1, drop = FALSE] else NULL
    }
    if (survMod == "spline-PH" && !is.null(W) && ncol(W) == 1)
        sims.list$gammas <- sims.list$gammas[, 1, drop = FALSE]
    colnames(sims.list$betas) <- colnames(X)
    colnames(sims.list$gammas) <- if (survMod == "weibull-PH") {
        if (is.null(W)) "(Intercept)" else c("(Intercept)", colnames(W)) 
    } else colnames(W)
    if (survMod == "spline-PH")
        colnames(sims.list$Bs.gammas) <- paste("gamma.bs", seq_len(ncW2), sep = "")
    colnames(sims.list$D) <- paste("D[", row(VC)[lower.tri(VC, TRUE)], ", ", 
        col(VC)[lower.tri(VC, TRUE)], "]", sep = "") 
    sims.list <- sims.list[!sapply(sims.list, is.null)]
    out$modes <- lapply(sims.list, function (x) {
        m <- function (x) {
            d <- density(x, bw = "nrd", adjust = 3, n = 1000)
            d$x[which.max(d$y)]
        }
        if (is.matrix(x)) apply(x, 2, m) else m(x)
    })
    out$coefficients <- lapply(sims.list, 
        function (x) if (is.matrix(x)) colMeans(x) else mean(x))
    out$StErr <- lapply(sims.list, 
        function (x) {
            f <- function (x) {
                acf.x <- drop(acf(x, lag.max = 0.4*length(x), plot = FALSE)$acf)[-1]
                acf.x <- acf.x[seq_len(rle(acf.x > 0)$lengths[1])]
                ess <- length(x) / (1 + 2 * sum(acf.x))
                sqrt(var(x) / ess)
            }
            if (is.matrix(x)) apply(x, 2, f) else f(x)
    })
    out$StDev <- lapply(sims.list, function (x) 
        if (is.matrix(x)) apply(x, 2, sd) else sd(x))
    out$CIs <- lapply(sims.list, 
        function (x) if (is.matrix(x)) apply(x, 2, quantile, probs = c(0.025, 0.975)) 
            else quantile(x, probs = c(0.025, 0.975)))
    out$vcov <- var(do.call(cbind, sims.list)) 
    D <- matrix(0, nrow(VC), ncol(VC))
    D[lower.tri(D, TRUE)] <- out$coefficients$D
    D <- D + t(D)
    diag(D) <- diag(D) / 2
    out$coefficients$D <- D
    dimnames(out$coefficients$D) <- dimnames(VC)
    D <- matrix(0, nrow(VC), ncol(VC))
    D[lower.tri(D, TRUE)] <- out$modes$D
    D <- D + t(D)
    diag(D) <- diag(D) / 2
    out$modes$D <- D
    dimnames(out$modes$D) <- dimnames(VC)    
    if (one.RE) {
        out$coefficients$D <- out$coefficients$D[1, 1, drop = FALSE]
        dimnames(out$coefficients$D) <- list(colnames(Z)[1], colnames(Z)[1])
        if (param == "shared-RE") {
            out$coefficients$alphas <- out$coefficients$alphas[1]
            names(out$coefficients$alphas) <- colnames(Z)[1]
        }
        x$Z <- x$Z[, 1, drop = FALSE]
        if (!is.null(x$Ztime))
            x$Ztime <- x$Ztime[, 1, drop = FALSE]
    }
    list.thetas <- out$coefficients
    thetas <- unlist(as.relistable(list.thetas))
    D.hat <- h(postMeans, thetas)
    M <- nrow(sims.list$betas)
    hat.Ds <- numeric(M)
    for (m in seq_len(M)) {
        postMeans.m <- matrix(ranef[m, ], ncol = ncol(Z), byrow = TRUE)
        thetas.m <- lapply(sims.list, function (x)
            if (is.matrix(x)) x[m, ] else x[m])
        D <- matrix(0, nrow(VC), ncol(VC))
        D[lower.tri(D, TRUE)] <- thetas.m$D
        D <- D + t(D)
        diag(D) <- diag(D) / 2
        thetas.m$D <- D
        thetas.m <- unlist(as.relistable(thetas.m))
        hat.Ds[m] <- h(postMeans.m, thetas.m)
    }
    hat.D <- mean(hat.Ds, na.rm = TRUE)
    out$pD <- hat.D - D.hat 
    out$DIC <- out$pD + hat.D
    rm(ranef, Bs, codaFit, sims.list)
    out$x <- x
    out$y <- y
    out$id <- id
    names(out$id) <- factor(lmeObject$groups[[1]])
    out$times <- data[[timeVar]]
    out$data <- data
    out$data.id <- data.id
    out$survMod <- survMod
    out$termsYx <- TermsX
    out$termsYz <- TermsZ
    if (param %in% c("td-extra", "td-both")) {
        out$termsYx.deriv <- TermsX.deriv
        out$termsYz.deriv <- TermsZ.deriv
    }
    out$termsT <- survObject$terms
    out$formYx <- formYx
    out$formYz <- formYz
    out$formT <- formT
    out$timeVar <- timeVar
    out$control <- con
    out$param <- param
    out$extraForm <- extraForm
    out$robust <- robust
    out$robust.b <- robust.b
    out$df <- df
    out$df.b <- df.b
    out$priors <- Data[grep("prior", names(Data), fixed = TRUE)]
    out$call <- cl
    class(out) <- "JMbayes"
    out
}
