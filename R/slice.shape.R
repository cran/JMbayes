slice.shape <-
function (logPost, current.shape, current.logPost, step, iters = 6) {
    logy <- current.logPost - rexp(1)
    tt <- runif(1, 0, step)
    L <- current.shape - tt
    R <- current.shape + step - tt
    while (L > 0 && logPost(L)[[1]] > logy) {
        L <- L - step
    }
    while (logPost(R)[[1]] > logy) {
        R <- R + step
    }
    L <- max(0, L)
    count <- 0
    repeat {
        count <- count + 1
        new.shape <- runif(1, L, R)
        new.logPost <- logPost(new.shape)[[1]]
        if (new.logPost >= logy || count > iters)
            break
        if (new.shape < current.shape) {
            L <- new.shape
        } else {
            R <- new.shape
        }
    }
    c(list(new.shape = new.shape, fail = count > iters), logPost(new.shape))
}
