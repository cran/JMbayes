# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("Dynamic Predictions using Joint Models"),
    
    sidebarPanel(
        wellPanel(
            fileInput('RDfile', 'Load the R Workspace with the fitted joint model',
                      accept = NULL),
            
            fileInput('patientFile', 'Load subject data',
                      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),

            withTags(div(class = 'row-fluid',
                         div(class = 'span3', radioButtons('sep', 'Separator', c(Comma = ',', 
                                                                                 Semicolon = ';', Tab = '\t'), ',')),
                         div(class = 'span3', radioButtons('dec', 'Decimal', c(Dot = '.', Comma = ','), '.')),
                         div(class = 'span5', radioButtons('quote', 'Quote', c(None = '', 'Double Quote' = '"', 
                                                                               'Single Quote' = "'"), '"'))
            ))
        ),
        
        wellPanel(
            uiOutput("modelChoose"),            
            
            uiOutput("obsChoose"),
            
            radioButtons('TypePlot', 'Type of Plot', 
                         c("Survival" = 'surv', "Cumulative Incidence" = "cumInc",
                           "Stick Man" = 'stickMan', "Longitudinal" = 'longitudinal'), 
                         "surv"),
            
            withTags(div(class = 'row-fluid',
                         div(class = 'span7', numericInput("time", "Target horizon time:", NULL)),
                         div(class = 'span5', checkboxInput("extra", "Add horizon time to the last visit time", FALSE))
            )),
            
            uiOutput("lastTime"),
            
            numericInput("M", "Monte Carlo samples:", 200),
            
            withTags(div(class = 'row-fluid',
                         div(class = 'span5', downloadButton('downloadData', 'Download Event-free Probabilities')),
                         div(class = 'span5', downloadButton('downloadPlot', 'Download Plot'))
            ))
        )
    ),
    
    mainPanel(
        tabsetPanel(
            tabPanel("Data", tableOutput('contents'), uiOutput("message")),
            tabPanel("Event-free Probabilities", tableOutput('survprobs'), uiOutput("message2")),
            tabPanel("Plot", plotOutput('plot'))
        )
    )
))
