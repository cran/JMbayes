replace.inprod <-
function (body.model, Data, param, extraForm, program) {
    mt <- deparse(body.model, width.cutoff = 200L)
    # long
    ncX <- Data$ncX
    rplc <- paste(paste("betas[", 1:ncX, "] * X[j, ", 1:ncX, "]", sep = ""), collapse = " + ")
    mt <- gsub("inprod(betas[1:ncX], X[j, 1:ncX])", rplc, mt, fixed = TRUE)
    ncZ <- Data$ncZ
    rplc <- paste(paste("b[i, ", 1:ncZ, "] * Z[j, ", 1:ncZ, "]", sep = ""), collapse = " + ")
    mt <- gsub("inprod(b[i, 1:ncZ], Z[j, 1:ncZ])", rplc, mt, fixed = TRUE)
    # surv
    ncW <- Data$ncW
    rplc <- paste(paste("gammas[", 1:ncW, "] * W[i, ", 1:ncW, "]", sep = ""), collapse = " + ")
    mt <- gsub("inprod(gammas[1:ncW], W[i, 1:ncW])", rplc, mt, fixed = TRUE)
    if (!is.null(Data$ncW2)) {
        ncW2 <- Data$ncW2
        rplc <- paste(paste("Bs.gammas[", 1:ncW2, "] * W2[i, ", 1:ncW2, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(Bs.gammas[1:ncW2], W2[i, 1:ncW2])", rplc, mt, fixed = TRUE)
        rplc <- paste(paste("Bs.gammas[", 1:ncW2, "] * W2s[K * (i - 1) + k, ", 1:ncW2, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(Bs.gammas[1:ncW2], W2s[K * (i - 1) + k, 1:ncW2])", rplc, mt, fixed = TRUE)        
    }
    # surv-value
    if (param %in% c("td-value", "td-both")) {
        rplc <- paste(paste("betas[", 1:ncX, "] * Xtime[i, ", 1:ncX, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(betas[1:ncX], Xtime[i, 1:ncX])", rplc, mt, fixed = TRUE)
        rplc <- paste(paste("b[i, ", 1:ncZ, "] * Ztime[i, ", 1:ncZ, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(b[i, 1:ncZ], Ztime[i, 1:ncZ])", rplc, mt, fixed = TRUE)
        #
        rplc <- paste(paste("betas[", 1:ncX, "] * Xs[K*(i - 1) + k, ", 1:ncX, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(betas[1:ncX], Xs[K * (i - 1) + k, 1:ncX])", rplc, mt, fixed = TRUE)
        rplc <- paste(paste("b[i, ", 1:ncZ, "] * Zs[K*(i - 1) + k, ", 1:ncZ, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(b[i, 1:ncZ], Zs[K * (i - 1) + k, 1:ncZ])", rplc, mt, fixed = TRUE)
    }
    if (param %in% c("td-extra", "td-both")) {
        indFixed <- extraForm$indFixed
        indRandom <- extraForm$indRandom
        ncX.deriv <- Data$ncX.deriv
        ncZ.deriv <- Data$ncZ.deriv
        rplc <- paste(paste("betas[", indFixed, "] * Xtime.deriv[i, ", 1:ncX.deriv, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(betas[indFixed], Xtime.deriv[i, 1:ncX.deriv])", rplc, mt, fixed = TRUE)
        rplc <- paste(paste("b[i, ", indRandom, "] * Ztime.deriv[i, ", 1:ncZ.deriv, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(b[i, indRandom], Ztime.deriv[i, 1:ncZ.deriv])", rplc, mt, fixed = TRUE)
        #
        rplc <- paste(paste("betas[", indFixed, "] * Xs.deriv[K * (i - 1) + k, ", 1:ncX.deriv, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(betas[indFixed], Xs.deriv[K * (i - 1) + k, 1:ncX.deriv])", rplc, mt, fixed = TRUE)
        rplc <- paste(paste("b[i, ", indRandom, "] * Zs.deriv[K * (i - 1) + k, ", 1:ncZ.deriv, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(b[i, indRandom], Zs.deriv[K * (i - 1) + k, 1:ncZ.deriv])", rplc, mt, fixed = TRUE)
    }
    if (param == "shared-RE") {
        nb <- Data$nb
        rplc <- paste(paste("alphas[", 1:nb, "] * b[i, ", 1:nb, "]", sep = ""), collapse = " + ")
        mt <- gsub("inprod(alphas[1:nb], b[i, 1:nb])", rplc, mt, fixed = TRUE)        
    }
    mt <- gsub("%_%", "", mt)
    if (program %in% c("OpenBUGS", "openbugs") && param == "shared-RE") {
        rplc <- "evTime[i] ~ dweib(sigma.t, mut[i]) C(censTime[i], )"
        mt <- gsub("evTime[i] ~ dweib(sigma.t, mut[i])  I(censTime[i], )", rplc, mt, fixed = TRUE)
    }
    if (program %in% c("JAGS", "jags") && param == "shared-RE") {
        rplc <- "is.censored[i] ~ dinterval(evTime[i], censTime[i])\n\tevTime[i] ~ dweib(sigma.t, mut[i])"
        mt <- gsub("evTime[i] ~ dweib(sigma.t, mut[i])  I(censTime[i], )", rplc, mt, fixed = TRUE)
    }
    c("model", mt)
}
