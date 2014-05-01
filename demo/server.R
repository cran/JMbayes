shinyServer(function(input, output) {    
    findJMbayesObjs <- reactive({
        inFile <- input$RDfile
        load(inFile$datapath)
        objs <- ls()
        ss <- sapply(objs, function (o) class(get(o))) == "JMbayes"
        if (all(!ss)) {
            stop("\nIt seems that there is no joint model fitted by jointModelBayes() ", 
                 "in the workspace that you loaded ...")
        }
        out <- objs[ss]
        names(out) <- out
        out
    })
    
    output$modelChoose <- renderUI({
        if (!is.null(input$RDfile)) {
            inFile <- input$RDfile
            load(inFile$datapath)
            objs <- findJMbayesObjs()
            selectInput("model", "Choose joint model:", objs, objs[1L])
        }
    })
    
    loadObject <- reactive({
        if (!is.null(input$RDfile) && !is.null(input$model)) {
            inFile <- input$RDfile
            load(inFile$datapath)
            objs <- findJMbayesObjs()
            choice <- input$model
            get(objs[input$model])
        }
    })
        
    dataObject <- reactive({
        if (!is.null(input$RDfile) && !is.null(input$model)) {
            object <- loadObject()
            forms <- object$Forms[c("formYx", "formYz", "formT")]
            forms$formT <- reformulate(attr(delete.response(object$Terms$termsT), "term.labels"))
            nams <- unique(unlist(lapply(forms, all.vars), use.names = FALSE))
            d <- object$Data$data[1:3, c(nams)]
            d[] <- lapply(d, function (x) {x[] <- NA; x})
            d
        }
    })
    
    ND <- reactive({
        if (!is.null(input$RDfile) && !is.null(input$model)) {
            if (is.null(input$patientFile)) {
                dataObject()
            } else {
                d <- dataObject()
                indF <- sapply(d, is.factor)
                levelsF <- lapply(d[indF], levels) 
                tt <- try(inData <- read.csv(input$patientFile$datapath, sep = input$sep, 
                                             quote=input$quote, colClasses = sapply(d, class),
                                             dec = input$dec), TRUE)
                if (!inherits(tt, "try-error")) {
                    inData <- cbind(inData, id = rep(1, nrow(inData)))
                    f <- function (f, l) {levels(f) <- l; f}
                    inData[indF] <- mapply(f, inData[indF], levelsF, SIMPLIFY = FALSE)
                    inData
                }
            }
        }
    })
    
    output$obsChoose <- renderUI({
        if (!is.null(input$patientFile)) {
            nr <- nrow(ND())
            if (!is.null(nr) && nr > 1) {
                sliderInput("obs", "Number of observations to use in prediction:", 
                            min = 1, max = nr, value = 1, step = 1,
                            animate = animationOptions(loop = TRUE))
            }
        }
    })
    
    output$lastTime <- renderUI({
        if (!is.null(input$patientFile)) {
            object <- loadObject()
            nd <- ND()
            times <- nd[[object$timeVar]]
            if (!is.null(times))
                numericInput("lasttime", "Last time point without event:", 
                             round(max(times, na.rm = TRUE), 2))
        }
    })
    
    output$message <- renderPrint({
        if (is.null(input$RDfile)) {
            cat("<br /><h3> <span style='color:black'> Welcome to the web interface",
                  "for producing dynamic predictions from joint models using package",
                  "<a href='http://cran.r-project.org/package=JMbayes'",
                  "style='text-decoration:none;' target='_blank'><b>JMbayes</b></a></h3>",
                "<br /><br /><h4>Load the R workspace containing the fitted joint model using the",
                "menus on the left to continue ...</h4>")
        } else {
            if (is.null(input$patientFile)) {
                dataset <- ND()
                classes <- sapply(dataset, class)
                classes[classes == "numeric"] <- "a numeric (continuous) variable."
                ind <- classes == "factor"
                classes[ind] <- paste("a factor (categorical) variable, with levels", 
                    sapply(dataset[ind], function (x) paste0("<i>", levels(x), "</i>", collapse = ", ")))
                msg1 <- "The data of the new subject should be stored in a CSV file with columns:"
                msg2 <- "<ul style='list-style-type:circle'>"
                msg3 <- paste("<li><b>", names(dataset), "</b>:", classes, "</li>", collapse = " ")
                msg4 <- "</ul>"
                msg5 <- "<br />You can use as a template the table above (e.g., copy-paste it in Excel and save it as CSV)." 
                msg6 <- "<u>Note:</u> R is case sensitive"
                cat("<br /><br /><p>", paste(msg1, msg2, msg3, msg4, msg5, msg6), "</p>")
            } else {
                d <- dataObject()
                tt <- try(inData <- read.csv(input$patientFile$datapath, sep = input$sep, 
                                             colClasses = sapply(d, class), dec = input$dec), TRUE)
                if (inherits(tt, "try-error"))
                    cat("<br /><span style='color:red'>Something went wrong when loading the data of the new subject;</span>",
                        "<br>try tweaking the 'Separator' and 'Decimal' options on the left panel ...")
            }
        } 
    })
    
    output$message2 <- renderPrint({
        if (is.null(input$RDfile) || is.null(input$patientFile))
            cat("<br />No subject data have been loaded yet ...")
    })
    
    sfits <- reactive({
        if (!is.null(input$patientFile) && input$TypePlot != 'longitudinal') {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            lastTimeUser <- input$lasttime
            lastTimeData <- max(nd[[object$timeVar]])
            sfits <- vector("list", n)
            for (i in 1:n) {
                lt <- if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData) 
                    lastTimeUser else NULL
                sfits[[i]] <- survfitJM(object, newdata = nd[1:i, ], M = input$M, 
                                        last.time = lt)
            }
            sfits
        }
    })
    
    lfits <- reactive({
        if (!is.null(input$patientFile) && input$TypePlot == 'longitudinal') {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            lastTimeUser <- input$lasttime
            lastTimeData <- max(nd[[object$timeVar]])
            lfits <- vector("list", n)
            for (i in 1:n) {
                lt <- if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData) 
                    lastTimeUser else NULL
                lfits[[i]] <- predict(object, newdata = nd[1:i, ], M = input$M, 
                                      returnData = TRUE, type = "Subject", 
                                      interval = "prediction", last.time = lt)
            }
            lfits
        }
    })
    
    sfits2 <- reactive({
        if (!is.na(input$time)) {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            lastTimeUser <- input$lasttime
            sfits <- vector("list", n)
            for (i in 1:n) {
                lastTimeData <- max(nd[1:i, object$timeVar])
                target.time <- if (input$extra) lastTimeData + input$time else input$time
                lt <- if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData) 
                    lastTimeUser else NULL
                if (i == n && !is.na(lastTimeUser) && lastTimeUser > lastTimeData && input$extra)
                    target.time <- lastTimeUser + input$time
                if (target.time > lastTimeData)
                    sfits[[i]] <- survfitJM(object, newdata = nd[1:i, ], 
                                            survTimes = target.time, M = input$M,
                                            last.time = lt)
            }
            sfits
        }
    })
    
    output$contents <- renderTable({
        if (!is.null(input$RDfile)) {
            nd <- ND()
            if (!is.null(input$patientFile)) {
                tt <- try(nd[, head(1:ncol(nd), -1)], TRUE)
                if (!inherits(tt, "try-error"))
                    nd <- tt
            }
            nd
        }
    })
    
    sprobs <- reactive({
        if (!is.null(input$patientFile)) {
            if (is.na(input$time)) {
                sfits. <- sfits()
            } else {
                sfits. <- sfits2()
            }
            nn <- if(is.na(input$obs)) length(sfits.) else input$obs
            ss <- sfits.[[nn]]
            if (!is.null(ss)) {
                sf <- sfits.[[nn]]
                f <- function (d, t) {
                    dd <- d[1, , drop = FALSE]
                    dd[1, ] <- c(as.vector(t), rep(1, ncol(dd) - 1))
                    round(rbind(dd, d), 4)
                }
                d <- mapply(f, sf$summaries, sf$last.time, SIMPLIFY = FALSE)
                d <- d[[1L]][, c(1,2,4,5), drop = FALSE]
                colnames(d) <- c("time", "Surv", "95% low", "95% upp")
                nr <- nrow(d)
                out <- if (nr > 10) d[round(seq(1, nr, length.out = 10)), ] else d
                rownames(out) <- NULL
                out
            } else NULL
        }
    })
    
    output$survprobs <- renderTable({
        sprobs()
    })
    
    output$plot <- renderPlot({
        if (!is.null(input$patientFile)) {
            if (input$TypePlot != 'longitudinal') {
                if (!is.na(input$time) && input$extra && input$TypePlot == "stickMan") {
                    sfits. <- sfits2()
                    nn <- if(is.na(input$obs)) length(sfits.) else input$obs
                    ss <- round(100 * (1 - sfits.[[nn]][["summaries"]][[1]][1, 2]))
                    draw.stick <- JMbayes:::draw.stick
                    xx <- seq(0.2, 1.1, 0.1)
                    yy <- seq(0.0, 1.8, 0.2)
                    cords <- expand.grid(xx, rev(yy))
                    cols <- rep("blue", 100)
                    cols[seq_len(ss)] <- "red"
                    op <- par(mar = c(0, 0, 2.1, 0))
                    plot(c(.25, 1.25), c(0, 2), type = "n", xaxt = 'n', yaxt = 'n', ann = FALSE)
                    title(paste0("Risk = ", ss, "%"), cex.main = 1.7)
                    for (cc in 1:100) {
                        draw.stick(cords[cc, 1], cords[cc, 2], linecol = cols[cc], scale = 0.2)
                    }
                    par(op)
                } else {
                    sfits. <- sfits()
                    nn <- if(is.na(input$obs)) length(sfits.) else input$obs
                    if (input$TypePlot == "surv") {
                        plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, fill.area = TRUE, 
                             include.y = TRUE, lwd = 2, ask = FALSE, cex = 2, main = "")                
                    }
                    if (input$TypePlot == "cumInc") {
                        plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, fill.area = TRUE, 
                             include.y = TRUE, lwd = 2, ask = FALSE, cex = 2, main = "",
                             fun = function (s) 1 - s, ylab = "Cumulative Incidence")                
                    }
                    if (input$extra) {
                        object <- loadObject()
                        nd <- ND()
                        nr <- nrow(nd)
                        lastTimeUser <- input$lasttime
                        lastTimeData <- max(nd[1:nn, object$timeVar])
                        target.time <- if (nn == nr && !is.na(lastTimeUser) && 
                                               lastTimeUser > lastTimeData) {
                            lastTimeUser + input$time
                        } else {
                            lastTimeData + input$time
                        }
                        abline(v = target.time, lty = 2, col = 2, lwd = 2)
                    }
                }
            } else {
                object <- loadObject()
                lfits. <- lfits()
                require("lattice")
                nn <- if(is.na(input$obs)) length(lfits.) else input$obs
                timeVar <- object$timeVar
                resp <- paste(object$Forms$formYx)[2L]
                tv <- lfits.[[nn]][[timeVar]]
                lastTimeData <- with(lfits.[[nn]], tv[!is.na(low)][1])
                lastTimeUser <- input$lasttime
                lastTime <- if (nn == length(lfits.)) 
                    max(lastTimeData, lastTimeUser, na.rm = TRUE) else lastTimeData
                form <- as.formula(paste("pred + low + upp +", resp, "~", timeVar))
                yl <- range(c(data.matrix(do.call(rbind, lfits.)[c("pred", "low", "upp")])),
                            na.rm = TRUE)
                yl <- yl + c(-0.05, 0.05) * yl
                target.time <- if (input$extra) {
                    object <- loadObject()
                    nd <- ND()
                    nr <- nrow(nd)
                    lastTimeUser <- input$lasttime
                    lastTimeData <- max(nd[1:nn, object$timeVar])
                    if (nn == nr && !is.na(lastTimeUser) && 
                                           lastTimeUser > lastTimeData) {
                        lastTimeUser + input$time
                    } else {
                        lastTimeData + input$time
                    }
                } else NA
                print(xyplot(form, data = lfits.[[nn]], last.time = lastTime, target.time = target.time,
                             panel = function (..., last.time, target.time) {
                                 xx <- ..1; yy <- ..2; ind <- ..4
                                 xx <- do.call(cbind, split(xx, ind))
                                 yy <- do.call(cbind, split(yy, ind))
                                 na.ind <- is.na(yy[, 2])
                                 lx <- xx[!na.ind, 2]; ll <- yy[!na.ind, 2]; lu <- yy[!na.ind, 3]
                                 lpolygon(c(lx, rev(lx)), c(ll, rev(lu)), border = "transparent", 
                                          col = "lightgrey")
                                 panel.xyplot(xx[, 1], yy[, 1], type = "l", lty = 1, col = 2, lwd = 2)
                                 panel.xyplot(xx[, 2], yy[, 2], type = "l", lty = 2, col = 1, lwd = 2)
                                 panel.xyplot(xx[, 3], yy[, 3], type = "l", lty = 2, col = 1, lwd = 2)
                                 panel.xyplot(xx[na.ind, 4], yy[na.ind, 4], type = "p", pch = 8, 
                                              cex = 1.1, col = 1)
                                 panel.abline(v = last.time, lty = 3)
                                 if (!is.na(target.time))
                                     panel.abline(v = target.time, lty = 2, col = 2, lwd = 2)
                             }, xlab = "Time", ylim = yl,
                             ylab = paste("Predicted", resp)))
            }
        }
    })
        
    output$downloadData <- downloadHandler(
        filename = function () { paste(input$dataset, '.csv', sep = '') },
        content = function (file) {
            write.table(sprobs(), file, sep = input$sep, dec = input$dec,
                        row.names = FALSE)
        }
    )
    
    output$downloadPlot <- downloadHandler(
        filename = function () { paste(input$dataset, '.png', sep = '') },
        content = function (file) {
            png(file)
            sfits. <- sfits()
            nn <- if(is.na(input$obs)) length(sfits.) else input$obs
            plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, fill.area = TRUE, 
                 include.y = TRUE, lwd = 2, ask = FALSE, cex = 2, main = "")
            dev.off()
        }
    )
})