tve <- function (x, df = NULL, knots = NULL, ord = 3) {
    if (is.null(df) && is.null(knots)) {
        stop("either 'df' or 'knots' need to be specified.\n")
    }
    if (is.null(knots) && !is.null(df)) {
        eps <- 0.1 * sd(x)
        min_x <- min(x) - eps
        max_x <- max(x) + eps
        dx <- (max_x - min_x) / df
        knots <- seq(min_x - (ord-1) * dx, max_x + (ord-1) * dx, by = dx)
    }
    out <- splines::splineDesign(knots, x, ord = ord, outer.ok = TRUE)
    attr(out, 'knots') <- knots
    attr(out, 'ord') <- ord
    attr(out, "class") <- c("tve", "basis", "matrix")
    out
}

makepredictcall.tve <- function (var, call) {
    if (as.character(call)[1L] != "tve") 
        return(call)
    at <- attributes(var)[c("knots", "ord")]
    xxx <- call[1L:2L]
    xxx[names(at)] <- at
    xxx
}