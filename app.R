library(shiny)

# Load the date once
# Load the mds data
mds_data <- data.frame(t(read.csv('data/mds_subspaces_pis_pk.csv', header = FALSE)))
subspace_names = c('pisani', 'pocket')
n_dims = paste0('dim', c(1:2))
# Name the columns dataframe
colnames(mds_data) <- paste(rep(subspace_names, each = length(n_dims)), n_dims, sep = '.')

# Plotting mds parameters
ax_mds <- list(linecolor = toRGB('black'), linewidth = 2, automargin = T, scaleratio =1, scaleanchor = "x",
               zerolinecolor = toRGB('grey'), zerolinewidth = 2, mirror = T, showline = T)
yax_mds <- ax_mds; yax_mds[['title']] <- 'Second Dimension'
xax_mds <- ax_mds; xax_mds[['title']] <- 'First Dimension'

ui <- fluidPage(

    # Application title
    titlePanel("CDK2: Consensus/ML Scorings"),
    
    fluidRow(wellPanel('Esta es una barra superior')),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fluidRow(
                column(6, selectInput(
                    inputId = 'mds_subspace',
                    label = 'cMDS subspace',
                    choices = c('Pisani Resiudes', 'Pocket Residues')
                ))
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                column(6, plotlyOutput(outputId = 'mdsPlot'))
                #column(6, plotOutput(outputId = 'distPlot_2'))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    #*** MDS Parameters ***
    mds_dataset <- reactive({
        mds_ <- switch (input$mds_subspace,
            'Pisani Resiudes' = mds_data[,1:2],
            'Pocket Residues' = mds_data[,3:4]
        )
        colnames(mds_) <- n_dims
        mds_
    })
    
    #***** OUTPUTS *****#
    output$mdsPlot <- renderPlotly({
        mds_ <- mds_dataset()
        mds_plot <- plot_ly(data = mds_, x = mds_$dim1, y = mds_$dim2)
        mds_plot <- mds_plot %>% layout(xaxis = xax_mds, yaxis = yax_mds, 
                                        dragmode = 'pan')
        mds_plot
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
