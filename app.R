library(shiny)

# Load the date once
# Load the RDS file
list_objs <- readRDS('data/data.rds')

mds_df <- list_objs$mds_df
# Plotting mds parameters
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
pal <- c("#f64b3c", "#5b8c5a", "#0f4c75", "#feb72b", "#29c7ac")
pal <- setNames(pal, c('active', 'inact_a', 'inact_b', 'inact_ope', 'dfg_out'))
pal_violin <- c('#2a7886', '#f1935c', '#ed7575')


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
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
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
                      color = mds_df[['Labels_conf']], 
                      text = text_hover,
                      size = mds_df$Inhib_mass,
                      hovertemplate = paste('%{text}')) %>% 
            layout(xaxis = xax_mds, yaxis = yax_mds, dragmode = 'pan')
        fig_mds
    })
    
    output$swarmPlot <- renderPlotly({
        # DkScore and DkLEff
        i = 0
        score_cols <- dk_score()
        fig_swarm <- plot_ly(type = 'violin', data = mds_df, 
                             boxpoints = "all", jitter = 1,
                             selectedpoints = c(1,2), text = text)
        for (score_ in score_cols) {
            fig_swarm <- add_trace(fig_swarm, y = mds_df[[score_]], pointpos = 0, 
                      marker = list(size = marker_size), 
                      name = strsplit(score_, '_')[[1]][1],
                      color = I(pal_violin[i + 1]))
        }
        fig_swarm <- layout(fig_swarm, xaxis = xax_sw, yaxis = yax_sw, dragmode = 'pan') %>% 
        highlight(on = "plotly_hover", off = "plotly_doubleclick")
        fig_swarm
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
