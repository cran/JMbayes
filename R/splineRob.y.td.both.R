splineRob.y.td.both <-
function () {
    for (i in 1:N) {
        # Longitudinal Part
        for (j in offset[i]:(offset[i+1] - 1)) {
            muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncZ], Z[j, 1:ncZ])
            y[j] ~ dt(muy[j], tau, df)
        }
        # Survival Part
        etaBaseline[i] <- inprod(gammas[1:ncW], W[i, 1:ncW])
        log.h0.T[i] <-  inprod(Bs.gammas[1:ncW2], W2[i, 1:ncW2])
        f.T[i] <- inprod(betas[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 1:ncZ], Ztime[i, 1:ncZ])
        f.T.deriv[i] <- inprod(betas[indFixed], Xtime.deriv[i, 1:ncX.deriv]) + inprod(b[i, indRandom], Ztime.deriv[i, 1:ncZ.deriv])
        log.hazard[i] <- log.h0.T[i] + etaBaseline[i] + alphas * f.T[i] + Dalphas * f.T.deriv[i]
        for (k in 1:K) {
            log.h0.s[i, k] <- inprod(Bs.gammas[1:ncW2], W2s[K * (i - 1) + k, 1:ncW2])
            f.s[i, k] <- inprod(betas[1:ncX], Xs[K*(i - 1) + k, 1:ncX]) + inprod(b[i, 1:ncZ], Zs[K*(i - 1) + k, 1:ncZ])
            f.s.deriv[i, k] <- inprod(betas[indFixed], Xs.deriv[K*(i - 1) + k, 1:ncX.deriv]) + inprod(b[i, indRandom], Zs.deriv[K*(i - 1) + k, 1:ncZ.deriv])
            SurvLong[i, k] <- wk[k] * exp(log.h0.s[i, k] + alphas * f.s[i, k] + Dalphas * f.s.deriv[i, k])
        }
        log.survival[i] <- - exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ])
        phi[i] <- C - (event[i] * log.hazard[i]) - log.survival[i]
        zeros[i] ~ dpois(phi[i])
        # Random Effects Part
        b[i, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
    }
    # Priors
    # Longitudinal Part
    betas[1:ncX] ~ dmnorm(priorMean.betas[], priorTau.betas[, ])
    tau ~ dgamma(priorA.tau, priorB.tau)
    # Survival Part
    gammas[1:ncW] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, ])
    alphas ~ dnorm(priorMean.alphas, priorTau.alphas)
    Dalphas ~ dnorm(priorMean.Dalphas, priorTau.Dalphas)
    Bs.gammas[1:ncW2] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, ])
    # Random Effects Part
    inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
