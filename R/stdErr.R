stdErr <-
function (x) {
    x <- as.matrix(x)
    vars <- apply(x, 2, var)
    ess <- effectiveSize(x)
    sqrt(vars / ess)
}
