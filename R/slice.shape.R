slice.shape <-
function (logPost, current.shapes, step, which, iters = 6L) {
    current.shape <- current.shapes[which]
    logy <- logPost(current.shape, which)[[1L]] - rexp(1L)
    tt <- runif(1L, 0, step)
    L <- current.shape - tt
    R <- current.shape + step - tt
    while (L > 0 && logPost(L, which)[[1L]] > logy) {
        L <- L - step
    }
    while (logPost(R, which)[[1L]] > logy) {
        R <- R + step
    }
    L <- max(0, L)
    count <- 0L
    repeat {
        count <- count + 1L
        new.shape <- runif(1L, L, R)
        new.logPost <- logPost(new.shape, which)
        if (new.logPost[[1L]] >= logy || count > iters)
            break
        if (new.shape < current.shape) {
            L <- new.shape
        } else {
            R <- new.shape
        }
    }
    c(list(new.shape = new.shape, fail = count > iters), new.logPost)
}
