library(shiny)
library(shinythemes)

# Load the date once
# Load the RDS file
list_objs <- readRDS('data/data.rds')

mds_df <- list_objs$mds_df
rownames(mds_df) <- c(1:dim(mds_df)[1])
# Plotting mds parameters
# Conf modebar
conf_ <- c(#"zoomIn2d", "zoomOut2d", #"select2d", 
          "toImage", "autoScale2d",
          "toggleSpikelines", "hoverCompareCartesian", 
          "hoverClosestCartesian") 
# Axes
ax_mds <- list(linecolor = toRGB('black'), automargin = T, 
               scaleratio =1, scaleanchor = "x",  dtick = 0.25, linewidth = 2,
               zerolinecolor = toRGB('grey'), zerolinewidth = 2, mirror = T, showline = T)
yax_mds <- ax_mds; yax_mds[['title']] <- 'Second Dimension'
xax_mds <- ax_mds; xax_mds[['title']] <- 'First Dimension'
# Violin axis
yax_sw <- ax_mds; yax_sw[['title']] <- 'ROC-AUC'; yax_sw[['range']] <- c(0.3, 1)
xax_sw <- ax_mds; xax_sw[['title']] <- 'Databases'
# Text hover
text_hover <- ~paste('<b>PDB id:</b>', mds_df$X,
                     '<br><b>Conf:</b>', mds_df$Labels_conf,
                     '<br><b>Ligand:</b>', mds_df$Inhib)
marker_size = 5
# Colors
pal <- c("#f64b3c", "#5b8c5a", "#2a7886", "#feb72b", "#29c7ac")
pal <- setNames(pal, c('active', 'inact_a', 'inact_b', 'inact_ope', 'dfg_out'))
pal_violin <- c('#2a7886', '#f1935c', '#ed7575')


ui <- fluidPage(
    theme=shinytheme('journal'),

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
                ))),
            fluidRow(
                column(6, selectInput(
                    inputId = 'dk_score',
                    label = 'Docking score type:',
                    choices = c('Docking Score', 'Ligand Efficiency Score')
                )))
        ),

        # Show a plot of the generated distribution
        mainPanel(
            fluidRow(
                column(6, plotlyOutput(outputId = 'mdsPlot')),
                column(6, plotlyOutput(outputId = 'swarmPlot'))
            ),
            verbatimTextOutput("click")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #*** MDS Parameters ***
    mds_subspaces <- reactive({
        switch (input$mds_subspace,
            'Pisani Resiudes' = c('pisani.dim1', 'pisani.dim2'),
            'Pocket Residues' = c('pocket.dim1', 'pocket.dim2')
        )
    })
    
    dk_score <- reactive({
        score_ <- switch (input$dk_score,
            'Docking Score' = 'dksc',
            'Ligand Efficiency Score' = 'dklef')
        dk_columns_ <- grep(score_, colnames(mds_df), ignore.case =T, value = T)
    })
   
    
    #####################
    #***** OUTPUTS *****#
    output$mdsPlot <- renderPlotly({
        mds_ <- mds_subspaces()
        fig_mds <- plot_ly(type = 'scatter', mode = 'markers', 
                           data = mds_df, colors = pal) %>% 
            add_trace(x = mds_df[[mds_[1]]], y = mds_df[[mds_[2]]],
                      key = rownames(mds_df),
                      color = mds_df[['Labels_conf']], 
                      text = text_hover,
                      size = mds_df$Inhib_mass,
                      hovertemplate = paste('%{text}')) %>% 
            add_trace(x = NA, y = NA, name = 'Selected') %>% # Ampty trace
            layout(xaxis = xax_mds, yaxis = yax_mds, dragmode = 'pan',
                   legend = list(title = list(text = '<b>Protein<br>Conformations:</b>')),
                   title = list(text = paste0('cMDS subspace (',
                                              input$mds_subspace, ')'), x = 0.1)) %>%
            config(modeBarButtonsToRemove = conf_, displaylogo = FALSE) %>%
            event_register("plotly_click")
        
    })
    
    observeEvent(event_data("plotly_click"), {
        d <- event_data("plotly_click")
        row_num <- d$key # Get the row number
        mds_ <- mds_subspaces()
        sliced_df <- mds_df[row_num, ]
        # Remove the previous selected point if it exist
        plotlyProxy('mdsPlot', session) %>%
            plotlyProxyInvoke("deleteTraces", -1)
        # Add the new trace
        plotlyProxy('mdsPlot', session) %>%
            plotlyProxyInvoke('addTraces', 
                  list(x = rep(sliced_df[[mds_[1]]], 2), 
                       y = rep(sliced_df[[mds_[2]]], 2),
                       size = 30, type = 'scatter', mode = 'markers',
                       marker = list(color = 'rgba(0, 0, 0, 0)',
                                     line = list(color = 'rgba(0, 0, 0, 1)',
                                     width = 2)), name = 'Selected'))
    }) 
    
    #**** Violin plots ****
    output$swarmPlot <- renderPlotly({
        # DkScore and DkLEff
        score_cols <- dk_score()

        fig_swarm <- plot_ly(type = 'violin', data = mds_df, 
                             points = "all", jitter = 1,
                             text = text)
        i = 0
        for (score_ in score_cols) {
            i = i + 1
            fig_swarm <- add_trace(fig_swarm, y = mds_df[[score_]], pointpos = 0, 
                      marker = list(size = marker_size), 
                      name = strsplit(score_, '_')[[1]][1],
                      color = I(pal_violin[i]),
                      box = list(visible = T), key = rownames(mds_df))
        }
        fig_swarm <- layout(fig_swarm, xaxis = xax_sw, yaxis = yax_sw, 
                            dragmode = 'select', 
                            legend = list(title = list(text = '<b>Ligand<br>Databases:</b>')),
                            title = list(text = paste0('AUC values (Vinardo ',
                                          input$dk_score, ')'), x = 0.1)) %>% 
            config(modeBarButtonsToRemove = conf_, displaylogo = FALSE) %>%
            event_register('plotly_selecting')
    })
    
    observeEvent(event_data("plotly_selected"), {
        d <- event_data("plotly_selected")
        row_num <- d$key # Get the row numbers
        mds_ <- mds_subspaces()
        sliced_df <- mds_df[row_num, ]
        # Remove the previous selected point if it exist
        plotlyProxy('mdsPlot', session) %>%
            plotlyProxyInvoke("deleteTraces", -1)
        # Add the new trace
        plotlyProxy('mdsPlot', session) %>%
            plotlyProxyInvoke('addTraces', 
                              list(x = sliced_df[[mds_[1]]], 
                                   y = sliced_df[[mds_[2]]],
                                   size = 30, type = 'scatter', mode = 'markers',
                                   marker = list(color = 'rgba(0, 0, 0, 0)',
                                                 line = list(color = 'rgba(0, 0, 0, 1)',
                                                             width = 2)), name = 'Selected'))
    }) 
    
    output$click <- renderPrint({
        d <- event_data("plotly_click")
        if (is.null(d)) "Click events appear here (double-click to clear)" else d
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
