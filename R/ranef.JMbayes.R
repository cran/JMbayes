ranef.JMbayes <-
function (object, postVar = FALSE, ...) {
    if (!inherits(object, "JMbayes"))
        stop("Use only with 'JMbayes' objects.\n")
    postMeans <- object$postMeans
    if (postVar) {
        attr(postMeans, "postVar") <- object$postVars
    }
    postMeans
}
