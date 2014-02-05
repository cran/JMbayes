shinyServer(function(input, output) {
    loadObject <- reactive({
        inFile <- input$RDfile
        load(inFile$datapath)
        objs <- ls()
        ss <- sapply(objs, function (o) class(get(o))) == "JMbayes"
        if (all(!ss))
            stop("\nIt seems that there is no joint model fitted by JMbayes in the workspace that you loaded ...")
        get(objs[ss][1])
    })
    
    dataObject <- reactive({
        if (!is.null(input$RDfile)) {
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
        if (!is.null(input$RDfile)) {
            if (is.null(input$patientFile)) {
                dataObject()
            } else {
                d <- dataObject()
                indF <- sapply(d, is.factor)
                levelsF <- lapply(d[indF], levels) 
                tt <- try(inData <- read.csv(input$patientFile$datapath, sep = input$sep, 
                                             quote=input$quote, colClasses = sapply(d, class)), TRUE)
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
                tt <- try(inData <- read.csv(input$patientFile$datapath, sep = input$sep, colClasses = sapply(d, class)), TRUE)
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
        if (!is.null(input$patientFile)) {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            sfits <- vector("list", n)
            for (i in 1:n) {
                sfits[[i]] <- survfitJM(object, newdata = nd[1:i, ], M = input$M)
            }
            sfits
        }
    })
    
    sfits2 <- reactive({
        if (!is.na(input$time)) {
            object <- loadObject()
            nd <- ND()
            n <- nrow(nd)
            sfits <- vector("list", n)
            for (i in 1:n) {
                if (input$time > max(nd[1:i, object$timeVar]))
                    sfits[[i]] <- survfitJM(object, newdata = nd[1:i, ], 
                                            survTimes = input$time, M = input$M)
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
            sfits. <- sfits()
            nn <- if(is.na(input$obs)) length(sfits.) else input$obs
            plot(sfits.[[nn]], estimator = "mean", conf.int = TRUE, fill.area = TRUE, 
                 include.y = TRUE, lwd = 2, ask = FALSE, cex = 2, main = "")
        }
    })
    
    output$downloadData <- downloadHandler(
        filename = function () { paste(input$dataset, '.csv', sep = '') },
        content = function (file) {
            write.csv(sprobs(), file)
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