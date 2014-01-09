fastSumID <-
function (x, group) {
    as.vector(rowsum.default(x, group, reorder = FALSE, na.rm = FALSE), mode = "numeric")
}
