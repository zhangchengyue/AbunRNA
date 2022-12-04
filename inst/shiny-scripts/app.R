library(shiny)


# load("abunMatrix.rda")
# load("countMatrix.rda")
# load("bigCond.rda")
# load("conditionsDF.rda")


# Define UI for random distribution app ----
ui <- fluidPage(

    # App title ----
    titlePanel("AbunRNA"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(


            # Input: Selector for choosing dataset ----
            selectInput(inputId = "dataset",
                        label = "Choose a dataset:",
                        choices = c("C.elegans twk-40 3 samples",
                                    "C.elegans twk-40 18 samples")),

            # Input: Numeric entry for number of obs to view ----
            numericInput(inputId = "obs",
                         label = "Number of observations to view:",
                         value = 10)

        ),

        mainPanel(

            tabsetPanel(type = "tabs",
                        tabPanel("Data Set", tableOutput("view")),
                        tabPanel("Heatmap", plotOutput("heatmap")),
                        tabPanel("PCA", tableOutput("pca")),
                        tabPanel("PCA Plot", plotOutput("pcaplot"))
            )

        )
    )
)



server <- function(input, output) {

    data(abunMatrix)
    data(countMatrix)
    data(bigCond)
    data(conditionsDF)

    datasetInput <- reactive({
        switch(input$dataset,
               "C.elegans twk-40 3 samples" = countMatrix,
               "C.elegans twk-40 18 samples" = abunMatrix)
    })

    output$view <- renderTable({
        head(datasetInput(), n = input$obs)
    })


    countHeat <- plotHeatMap(matrix = countMatrix, head = T)

    abunHeat <- plotHeatMap(matrix = abunMatrix, head = T)

    heatInput <- reactive({
        switch(input$dataset,
               "C.elegans twk-40 3 samples" = countHeat,
               "C.elegans twk-40 18 samples" = abunHeat)
    })


    output$heatmap <- renderPlot({
        heatInput()
    })


    countPlot <- plotPCA(matrix = countMatrix, scaleIt = TRUE,
                         conditions = conditionsDF,
                         col = "genotype")

    abunPlot <- plotPCA(matrix = abunMatrix, scaleIt = TRUE,
                        conditions = bigCond,
                        col = "genotype")

    pcaInput <- reactive({
        switch(input$dataset,
               "C.elegans twk-40 3 samples" = countPlot,
               "C.elegans twk-40 18 samples" = abunPlot)
    })


    output$pcaplot <- renderPlot({
        pcaInput()$Plot
    })

    output$pca <- renderTable({
        head(pcaInput()$PCA, n = input$obs)
    })

}

shinyApp(ui = ui, server = server)
