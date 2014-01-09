effectiveSize <-
function (x) {
    # copied (and made a bit more efficient) from the coda package
    spectrum0.ar <- function (x) {
        v0 <- numeric(ncol(x))
        z <- 1:nrow(x)
        res <- as.matrix(lm.fit(cbind(1, z), cbind(x, x))$residuals)
        for (i in 1:ncol(x)) {
            if (identical(all.equal(sd(res[, i]), 0), TRUE)) {
                v0[i] <- 0
            } else {
                ar.out <- ar(x[, i], aic = TRUE)
                v0[i] <- ar.out$var.pred / (1 - sum(ar.out$ar))^2
            }
        }
        v0
    }
    x <- as.matrix(x)
    spec <- spectrum0.ar(x)
    ifelse(spec == 0, 0, nrow(x) * apply(x, 2, var) / spec)
}
