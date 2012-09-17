write.model.JMbayes <-
function (model, con = "model.bug", Data, param, extraForm, program) {
    model <- replace.inprod(body(model), Data, param, extraForm, program)
    writeLines(model, con)
}
