splineRob.yb.sharedRE <-
function () {
    for (i in 1:N) {
        # Longitudinal Part
        for (j in offset[i]:(offset[i+1] - 1)) {
            muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncZ], Z[j, 1:ncZ])
            y[j] ~ dt(muy[j], tau, df)
        }
        # Survival Part
        log.h0.T[i] <-  inprod(Bs.gammas[1:ncW2], W2[i, 1:ncW2])
        etaBaseline[i] <- inprod(gammas[1:ncW], W[i, 1:ncW])
        f.T[i] <- inprod(alphas[1:nb], b[i, 1:nb])
        log.hazard[i] <- log.h0.T[i] + etaBaseline[i] + f.T[i]
        for (k in 1:K) {
            log.h0.s[i, k] <- inprod(Bs.gammas[1:ncW2], W2s[K * (i - 1) + k, 1:ncW2])
            SurvLong[i, k] <- wk[k] * exp(log.h0.s[i, k])
        }
        log.survival[i] <- - exp(etaBaseline[i] + f.T[i]) * P[i] * sum(SurvLong[i, ])
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
    alphas[1:nb] ~ dmnorm(priorMean.alphas[], priorTau.alphas[, ])
    Bs.gammas[1:ncW2] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, ])
    # Random Effects Part
    inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
