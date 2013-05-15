namesJMbayes <-
function (survMod, param, robust, robust.b) {
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
    namesSplRob.y.td.value <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSplRob.y.td.extra <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    #
    namesSplRob.y.td.both <- c(
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
    namesSplRob.y.sharedRE <- c(
    "N", "K", "df",
    "offset", "y", "X", "Z", "ncX",
    "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk", 
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Bs.gammas", "priorTau.Bs.gammas",
    "priorR.D", "priorK.D")
    namesSplRob.b.td.value <- c(
        "N", "K", "df.b",
        "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
        "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
        "priorR.D", "priorK.D")
    #
    namesSplRob.b.td.extra <- c(
        "N", "K", "df.b",
        "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
        "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
        "priorR.D", "priorK.D")
    #
    namesSplRob.b.td.both <- c(
        "N", "K", "df.b",
        "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
        "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
        "priorR.D", "priorK.D")
    #
    namesSplRob.b.sharedRE <- c(
        "N", "K", "df.b",
        "offset", "y", "X", "Z", "ncX",
        "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk", 
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Bs.gammas", "priorTau.Bs.gammas",
        "priorR.D", "priorK.D")
    #
    namesSplRob.yb.td.value <- c(
        "N", "K", "df", "df.b",
        "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
        "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
        "priorR.D", "priorK.D")
    #
    namesSplRob.yb.td.extra <- c(
        "N", "K", "df", "df.b",
        "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
        "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
        "priorR.D", "priorK.D")
    #
    namesSplRob.yb.td.both <- c(
        "N", "K", "df", "df.b",
        "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
        "event", "zeros", "W", "ncW", "W2", "W2s", "ncW2", "C", "P", "wk",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorMean.Bs.gammas", "priorTau.Bs.gammas",
        "priorR.D", "priorK.D")
    #
    namesSplRob.yb.sharedRE <- c(
        "N", "K", "df", "df.b",
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
    namesWeibRob.y.td.value <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeibRob.y.td.extra <- c(
    "N", "K", "df",
    "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
    "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeibRob.y.td.both <- c(
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
    namesWeibRob.y.sharedRE <- c(
    "N", "df",
    "offset", "y", "X", "Z", "ncX",
    "evTime", "censTime", "W", "ncW",
    "nb", "mu0",
    "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
    "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
    "priorR.D", "priorK.D")
    #
    namesWeibRob.b.td.value <- c(
        "N", "K", "df.b",
        "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
        "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    #
    namesWeibRob.b.td.extra <- c(
        "N", "K", "df.b",
        "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
        "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    #
    namesWeibRob.b.td.both <- c(
        "N", "K", "df.b",
        "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
        "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    #
    namesWeibRob.b.sharedRE <- c(
        "N", "df.b",
        "offset", "y", "X", "Z", "ncX",
        "evTime", "censTime", "W", "ncW",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    
    namesWeibRob.yb.td.value <- c(
        "N", "K", "df", "df.b",
        "offset", "y", "X", "Xtime", "Xs", "Z", "Ztime", "Zs", "ncX",
        "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    #
    namesWeibRob.yb.td.extra <- c(
        "N", "K", "df", "df.b",
        "offset", "y", "X", "Xtime.deriv", "Xs.deriv", "Z", "Ztime.deriv", "Zs.deriv", "ncX",
        "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    #
    namesWeibRob.yb.td.both <- c(
        "N", "K", "df", "df.b",
        "offset", "y", "X", "Xtime", "Xtime.deriv", "Xs", "Xs.deriv", "Z", "Ztime", 
        "Ztime.deriv", "Zs", "Zs.deriv", "ncX",
        "event", "logTime", "zeros", "W", "ncW", "C", "P", "wk", "st",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", 
        "priorMean.Dalphas", "priorTau.Dalphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    #
    namesWeibRob.yb.sharedRE <- c(
        "N", "df", "df.b",
        "offset", "y", "X", "Z", "ncX",
        "evTime", "censTime", "W", "ncW",
        "nb", "mu0",
        "priorMean.betas", "priorTau.betas", "priorA.tau", "priorB.tau",
        "priorMean.gammas", "priorTau.gammas", "priorMean.alphas", "priorTau.alphas", "priorA.sigma.t", "priorB.sigma.t",
        "priorR.D", "priorK.D")
    ###############################################
    switch(paste(survMod, param, robust, robust.b, sep = "/"),
        "weibull-PH/td-value/FALSE/FALSE" = namesWeib.td.value,
        "weibull-PH/td-extra/FALSE/FALSE" = namesWeib.td.extra,
        "weibull-PH/td-both/FALSE/FALSE" = namesWeib.td.both,
        "weibull-PH/shared-RE/FALSE/FALSE" = namesWeib.sharedRE,
        "spline-PH/td-value/FALSE/FALSE" = namesSpl.td.value,
        "spline-PH/td-extra/FALSE/FALSE" = namesSpl.td.extra,
        "spline-PH/td-both/FALSE/FALSE" = namesSpl.td.both,
        "spline-PH/shared-RE/FALSE/FALSE" = namesSpl.sharedRE,
        ### 
        "weibull-PH/td-value/TRUE/FALSE" = namesWeibRob.y.td.value,
        "weibull-PH/td-extra/TRUE/FALSE" = namesWeibRob.y.td.extra,
        "weibull-PH/td-both/TRUE/FALSE" = namesWeibRob.y.td.both,
        "weibull-PH/shared-RE/TRUE/FALSE" = namesWeibRob.y.sharedRE,
        "spline-PH/td-value/TRUE/FALSE" = namesSplRob.y.td.value,
        "spline-PH/td-extra/TRUE/FALSE" = namesSplRob.y.td.extra,
        "spline-PH/td-both/TRUE/FALSE" = namesSplRob.y.td.both,
        "spline-PH/shared-RE/TRUE/FALSE" = namesSplRob.y.sharedRE,
        ### 
        "weibull-PH/td-value/FALSE/TRUE" = namesWeibRob.b.td.value,
        "weibull-PH/td-extra/FALSE/TRUE" = namesWeibRob.b.td.extra,
        "weibull-PH/td-both/FALSE/TRUE" = namesWeibRob.b.td.both,
        "weibull-PH/shared-RE/FALSE/TRUE" = namesWeibRob.b.sharedRE,
        "spline-PH/td-value/FALSE/TRUE" = namesSplRob.b.td.value,
        "spline-PH/td-extra/FALSE/TRUE" = namesSplRob.b.td.extra,
        "spline-PH/td-both/FALSE/TRUE" = namesSplRob.b.td.both,
        "spline-PH/shared-RE/FALSE/TRUE" = namesSplRob.b.sharedRE,
        ###
        "weibull-PH/td-value/TRUE/TRUE" = namesWeibRob.yb.td.value,
        "weibull-PH/td-extra/TRUE/TRUE" = namesWeibRob.yb.td.extra,
        "weibull-PH/td-both/TRUE/TRUE" = namesWeibRob.yb.td.both,
        "weibull-PH/shared-RE/TRUE/TRUE" = namesWeibRob.yb.sharedRE,
        "spline-PH/td-value/TRUE/TRUE" = namesSplRob.yb.td.value,
        "spline-PH/td-extra/TRUE/TRUE" = namesSplRob.yb.td.extra,
        "spline-PH/td-both/TRUE/TRUE" = namesSplRob.yb.td.both,
        "spline-PH/shared-RE/TRUE/TRUE" = namesSplRob.yb.sharedRE)
}
