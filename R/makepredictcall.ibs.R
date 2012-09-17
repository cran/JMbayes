makepredictcall.ibs <-
function (var, call) {
    if (as.character(call)[1L] != "ibs") 
        return(call)
    at <- attributes(var)[c("knots", "Boundary.knots", "intercept")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}
