WeibRob.yb.sharedRE <-
function () {
    for (i in 1:N) {
        # Longitudinal Part
        for (j in offset[i]:(offset[i+1] - 1)) {
            muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncZ], Z[j, 1:ncZ])
            y[j] ~ dt(muy[j], tau, df)
        }
        # Survival Part
        log(mut[i]) <- inprod(gammas[1:ncW], W[i, 1:ncW]) + inprod(alphas[1:nb], b[i, 1:nb])
        evTime[i] ~ dweib(sigma.t, mut[i]) %_% I(censTime[i], )
        # Random Effects Part
        b[i, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
    }
    # Priors
    # Longitudinal Part
    betas[1:ncX] ~ dmt(priorMean.betas[], priorTau.betas[, ], df.b)
    tau ~ dgamma(priorA.tau, priorB.tau)
    # Survival Part
    gammas[1:ncW] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, ])
    alphas[1:nb] ~ dmnorm(priorMean.alphas[], priorTau.alphas[, ])
    sigma.t ~ dgamma(priorA.sigma.t, priorB.sigma.t)
    # Random Effects Part
    inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
