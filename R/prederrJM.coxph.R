prederrJM.coxph <-
function (object, newdata, Tstart, Thoriz, lossFun = c("absolute", "square"), 
                               interval = FALSE, idVar = "id", timeVar = "time", ...) {
    if (!inherits(object, "coxph"))
        stop("Use only with 'coxph' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(newdata[[idVar]]))
        stop("'idVar' not in 'newdata.\n'")
    lossFun <- if (is.function(lossFun)) {
        lf <- lossFun
        match.fun(lossFun)
    } else {
        lf <- match.arg(lossFun)
        if (lf == "absolute") function (x) abs(x) else function (x) x*x
    }
    id <- newdata[[idVar]]
    id <- match(id, unique(id))
    TermsT <- object$terms
    SurvT <- model.response(model.frame(TermsT, newdata)) 
    Time <- SurvT[, 1]
    newdata2 <- newdata[Time > Tstart, ]
    SurvT <- model.response(model.frame(TermsT, newdata2)) 
    Time <- SurvT[, 1]
    delta <- SurvT[, 2]
    aliveThoriz <- newdata2[Time > Thoriz & newdata2[[timeVar]] < Tstart, ]
    deadThoriz <- newdata2[Time <= Thoriz & delta == 1 & newdata2[[timeVar]] < Tstart, ]
    indCens <- Time < Thoriz & delta == 0 & newdata2[[timeVar]] < Tstart
    censThoriz <- newdata2[indCens, ]
    nr <- length(unique(newdata2[[idVar]])) 
    idalive <- unique(aliveThoriz[[idVar]])
    iddead <- unique(deadThoriz[[idVar]])
    idcens <- unique(censThoriz[[idVar]])
    idalive. <- factor(aliveThoriz[[idVar]], idalive)
    iddead. <- factor(deadThoriz[[idVar]], iddead)
    idcens. <- factor(censThoriz[[idVar]], idcens)
    aliveThoriz.id <- aliveThoriz[tapply(row.names(aliveThoriz), idalive., tail, 1), ]
    Surv.aliveThoriz <- c(summary(survfit(object, newdata = aliveThoriz.id), times = Thoriz)$surv)
    deadThoriz.id <- deadThoriz[tapply(row.names(deadThoriz), iddead., tail, 1), ]
    Surv.deadThoriz <- c(summary(survfit(object, newdata = deadThoriz.id), times = Thoriz)$surv)
    if (nrow(censThoriz)) {
        censThoriz.id <- censThoriz[tapply(row.names(censThoriz), idcens., tail, 1), ]
        Surv.censThoriz <- c(summary(survfit(object, newdata = censThoriz.id), times = Thoriz)$surv)
        tt <- model.response(model.frame(TermsT, censThoriz.id))[, 1]
        nn <- length(tt)
        weights <- numeric(nn)
        for (i in seq_len(nn)) {
            weights[i] <- c(summary(survfit(object, newdata = censThoriz.id[i, ]), times = Thoriz)$surv) /
                c(summary(survfit(object, newdata = censThoriz.id[i, ]), times = tt[i])$surv)
        }
    } else {
        Surv.censThoriz <- weights <- NA
    }
    prederr <- if (!interval) {
        (1/nr) * sum(lossFun(1 - Surv.aliveThoriz), lossFun(0 - Surv.deadThoriz),
                     weights * lossFun(1 - Surv.censThoriz) + (1 - weights) * lossFun(0 - Surv.censThoriz))
    } else {
        TimeCens <- model.response(model.frame(TermsT, newdata))[, 1]
        deltaCens <- 1 - model.response(model.frame(TermsT, newdata))[, 2]
        KMcens <- survfit(Surv(TimeCens, deltaCens) ~ 1)
        times <- TimeCens[TimeCens > Tstart & TimeCens < Thoriz & !deltaCens]
        times <- sort(unique(times))
        k <- as.numeric(table(times))
        w <- summary(KMcens, times = Tstart)$surv / summary(KMcens, times = times)$surv
        prederr.times <- sapply(times, 
                                function (t) prederrJM(object, newdata, Tstart, t,
                                                       interval = FALSE, idVar = idVar, simulate = simulate)$prederr)
        num <- sum(prederr.times * w * k, na.rm = TRUE)
        den <- sum(w * k, na.rm = TRUE)
        num / den
    }
    out <- list(prederr = prederr, nr = nr, Tstart = Tstart, Thoriz = Thoriz, interval = interval,
                classObject = class(object), nameObject = deparse(substitute(object)), lossFun = lf)
    class(out) <- "prederrJM"
    out
}
