namesJMbayes <-
function (survMod, param, robust) {
    namesSpl.td.value <- c(
    "N", "K",
    "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSpl.td.extra <- c(
    "N", "K",
    "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSpl.td.both <- c(
    "N", "K",
    "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSpl.sharedRE <- c(
    "N", "K",
    "offset", "y", "X", "Z", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk", 
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSplRob.td.value <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSplRob.td.extra <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSplRob.td.both <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSplRob.sharedRE <- c(
    "N", "K", "df",
    "offset", "y", "X", "Z", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk", 
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesWeib.td.value <- c(
    "N", "K",
    "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeib.td.extra <- c(
    "N", "K",
    "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeib.td.both <- c(
    "N", "K",
    "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeib.sharedRE <- c(
    "N",
    "offset", "y", "X", "Z", "ncX",
    "evTime", "censTime", "W", "ncW",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeibRob.td.value <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeibRob.td.extra <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeibRob.td.both <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeibRob.sharedRE <- c(
    "N", "df",
    "offset", "y", "X", "Z", "ncX",
    "evTime", "censTime", "W", "ncW",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    ###############################################
    switch(paste(survMod, param, robust, sep = "/"),
        "weibull-PH/td-value/FALSE" = namesWeib.td.value,
        "weibull-PH/td-extra/FALSE" = namesWeib.td.extra,
        "weibull-PH/td-both/FALSE" = namesWeib.td.both,
        "weibull-PH/shared-RE/FALSE" = namesWeib.sharedRE,
        "spline-PH/td-value/FALSE" = namesSpl.td.value,
        "spline-PH/td-extra/FALSE" = namesSpl.td.extra,
        "spline-PH/td-both/FALSE" = namesSpl.td.both,
        "spline-PH/shared-RE/FALSE" = namesSpl.sharedRE,
        "weibull-PH/td-value/TRUE" = namesWeibRob.td.value,
        "weibull-PH/td-extra/TRUE" = namesWeibRob.td.extra,
        "weibull-PH/td-both/TRUE" = namesWeibRob.td.both,
        "weibull-PH/shared-RE/TRUE" = namesWeibRob.sharedRE,
        "spline-PH/td-value/TRUE" = namesSplRob.td.value,
        "spline-PH/td-extra/TRUE" = namesSplRob.td.extra,
        "spline-PH/td-both/TRUE" = namesSplRob.td.both,
        "spline-PH/shared-RE/TRUE" = namesSplRob.sharedRE)
}
