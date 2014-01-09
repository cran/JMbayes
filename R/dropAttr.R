dropAttr <-
function (x) {
    if (is.null(x))
        return(NULL)
    if (is.matrix(x)) {
        d <- dim(x)
        x <- as.vector(x)
        dim(x) <- d
        x
    } else as.vector(x)
}
