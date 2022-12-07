library(shiny)
library(shinydashboard)
library(shinyjs)
library(ggfortify)
library(dplyr)


ui <- fluidPage(

    # App title ----
    titlePanel("AbunRNA"),

    sidebarLayout(
        sidebarPanel(
            useShinyjs(),
            box(id = "generate", width = "800px",
                title = "Demo or Run?",
                selectInput(inputId = "generate",
                            label = "Would you like to play with our demo or
                            generate your own matrix by uploading quantification
                            files?",
                            choices = c("Demo", "Upload files")
                )
            ),

            # Input: Selector for choosing dataset ----
            box(id = "demo", width = "800px",
                selectInput(inputId = "dataset",
                            label = "Choose a dataset:",
                            choices = c("C.elegans twk-40 3 samples",
                                        "C.elegans twk-40 18 samples"))
            ),

            box(id = "ownMatrix", width = "800px",
                fileInput(
                    inputId = "files",
                    label = "Choose Quantification Files",
                    multiple = TRUE),

                textInput(inputId = "species",
                          label = "Type in reference transcriptom organism.
                      e.g. Caenorhabditis elegans", value = "", width = NULL,
                          placeholder = NULL),

                box(id = "check",
                    width = "800px",
                    title = "Release Version",
                    selectInput(inputId = "check",
                                label = "Do you have specific version preference
                                on the reference transcriptom?",
                                choices = c(FALSE, TRUE))
                ),

                box(id = "release",
                    width = 80,
                    textInput(inputId = "release",
                              label = "Type in the release version below.
                          e.g. 107",
                              value = "")
                ),
                actionButton("go", "Run")
            ),

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


    data(countMatrix)
    data(abunMatrix)
    data(bigCond)
    data(conditionsDF)



    observeEvent(input$generate, {

        if (input$generate == "Demo") {
            shinyjs::hide(id = "ownMatrix")
            shinyjs::show(id = "demo")

            # countMatrix <- AbunRNA::countMatrix
            # abunMatrix <- AbunRNA::abunMatrix
            # bigCond <- AbunRNA::bigCond
            # conditionsDF <- AbunRNA::conditionsDF


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

        else {
            shinyjs::show(id = "ownMatrix")
            shinyjs::hide(id = "demo")

            output$view <- renderTable(NULL)
            output$heatmap <- renderPlot(NULL)
            output$pca <- renderTable(NULL)
            output$pcaplot <- renderTable(NULL)

            matrix <- data.frame()


            observeEvent(input$files, {
                print(input$files %>% dplyr::arrange(desc(size)))
            })

            observeEvent(input$check, {

                if (input$check == TRUE) {
                    shinyjs::show(id = "release")
                }
                else {
                    shinyjs::hide(id = "release")
                }
            })

            getMatrix <- eventReactive(input$go, {

                upload <- c()

                for(nr in seq_along(input$files)) {
                    upload[nr] <- input$files[[nr, 'datapath']]
                }

                matrix <- data.frame()
                if (length(upload) != 0) {

                    matrix <- generateMatrix(sfSeq = upload,
                                             species = input$species,
                                             release = input$release)

                } else {
                    ;
                }

                return(matrix)
            })

            getHeatmap <- eventReactive(input$go, {

                if (length(getMatrix()) != 0) {
                    heatMap <- plotHeatMap(matrix = getMatrix(), head = T)
                }
                else {
                    ;
                }
                return(heatmap)
            })


            getPCA <- eventReactive(input$go, {

                if (length(getMatrix()) != 0) {
                    pcaResult <- plotPCA(matrix = getMatrix(), scaleIt = TRUE)
                }
                else {
                    ;
                }
                return(pcaResult)
            })
            output$view <- renderTable({
                head(getMatrix(), n = input$obs)
                # head(analysis())
            })

            output$heatmap <- renderPlot({
                getHeatmap()
            })

            output$pcaplot <- renderPlot({
                getPCA()$Plot
            })

            output$pca <- renderTable({
                head(getPCA()$PCA, n = input$obs)
            })

        }
    })

}

shiny::shinyApp(ui, server)

# [END]
