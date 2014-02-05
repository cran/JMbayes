slice.shape <-
function (logPost, current.shape, step, iters = 6L) {
    logy <- logPost(current.shape)[[1L]] - rexp(1L)
    tt <- runif(1L, 0, step)
    L <- current.shape - tt
    R <- current.shape + step - tt
    while (L > 0 && logPost(L)[[1L]] > logy) {
        L <- L - step
    }
    while (logPost(R)[[1L]] > logy) {
        R <- R + step
    }
    L <- max(0, L)
    count <- 0L
    repeat {
        count <- count + 1L
        new.shape <- runif(1L, L, R)
        new.logPost <- logPost(new.shape)
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
