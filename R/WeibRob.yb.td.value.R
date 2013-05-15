WeibRob.yb.td.value <-
function () {
    for (i in 1:N) {
        # Longitudinal Part
        for (j in offset[i]:(offset[i+1] - 1)) {
            muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncZ], Z[j, 1:ncZ])
            y[j] ~ dt(muy[j], tau, df)
        }
        # Survival Part
        etaBaseline[i] <- inprod(gammas[1:ncW], W[i, 1:ncW])
        f.T[i] <- inprod(betas[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 1:ncZ], Ztime[i, 1:ncZ])
        log.hazard[i] <- log(sigma.t) + (sigma.t - 1) * logTime[i] + etaBaseline[i] + alphas * f.T[i]
        for (k in 1:K) {
            f.s[i, k] <- inprod(betas[1:ncX], Xs[K*(i - 1) + k, 1:ncX]) + inprod(b[i, 1:ncZ], Zs[K*(i - 1) + k, 1:ncZ])
            SurvLong[i, k] <- wk[k] * sigma.t * pow(st[i, k], sigma.t - 1) * exp(alphas * f.s[i, k])
        }
        log.survival[i] <- - exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ])
        phi[i] <- C - (event[i] * log.hazard[i]) - log.survival[i]
        zeros[i] ~ dpois(phi[i])
        # Random Effects Part
        b[i, 1:nb] ~ dmt(mu0[], inv.D[, ], df.b)
    }
    # Priors
    # Longitudinal Part
    betas[1:ncX] ~ dmnorm(priorMean.betas[], priorTau.betas[, ])
    tau ~ dgamma(priorA.tau, priorB.tau)
    # Survival Part
    gammas[1:ncW] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, ])
    alphas ~ dnorm(priorMean.alphas, priorTau.alphas)
    sigma.t ~ dgamma(priorA.sigma.t, priorB.sigma.t)
    # Random Effects Part
    inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
