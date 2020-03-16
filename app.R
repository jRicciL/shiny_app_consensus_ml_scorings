library(shiny)
library(shinythemes)
library(plotly)
library(dplyr)
library(shinyjs)

# Load the date once
df_consensus <- read.csv('data/cosensus_scorings.csv')
# Save it on a list to simulate the future behaviour
# scores = list(kmeans = df_consensus)
# Get the number of observations
n_ <- length(data) - 5 # Because always the first columns are not auc values
n_steps <- c(1: n_)

# Load the RDS file
list_objs <- readRDS('data/data.rds')

# Reference values
ref_values_dk_scores <- list(
    'dksc' = list('csar' = 0.8486, 'dud' = 0.675, 'dekois' = 0.7904),
    'dklef' = list('csar' = 0.9287, 'dud' = 0.7788, 'dekois' = 0.7342))

# funtion to use apply
get_col_names <- function(x, fig){
    trace_name <- paste(c(toupper(x[4]), x[5]), collapse = ' ')
}

#*********** MDS and Violin Plots ***************
mds_df <- list_objs$mds_df
rownames(mds_df) <- c(1:dim(mds_df)[1])
# Plotting mds parameters
# Conf modebar
conf_ <- c(#"zoomIn2d", "zoomOut2d", #"select2d", 
          "toImage", "autoScale2d",
          "toggleSpikelines", "hoverCompareCartesian", 
          "hoverClosestCartesian") 

#### LINE PLOT
ax_lp <- list(linecolor = toRGB('black'), mirror = T,
              linewidth = 3, showline = T, dtick = 25, marging = list(pad = 0),
              title = 'Number of used Protein Conformations')
yax_lp <- ax_lp; yax_lp[['range']] <- c(0.3, 1)
yax_lp[['dtick']] <- 0.1
yax_lp[['title']] <- 'ROC-AUC'
#### MDS and VIOLIN PLOT
# Axes
ax_mds <- list(linecolor = toRGB('black'), automargin = T, 
               scaleratio =1, scaleanchor = "x",  dtick = 0.25, linewidth = 3,
               zerolinecolor = toRGB('grey'), zerolinewidth = 2, mirror = T, showline = T)
yax_mds <- ax_mds; yax_mds[['title']] <- 'Second Dimension'
xax_mds <- ax_mds; xax_mds[['title']] <- 'First Dimension'
# Violin axis
yax_sw <- ax_mds; yax_sw[['title']] <- 'ROC-AUC'
yax_sw[['range']] <- c(0.3, 1)
yax_sw[['dtick']] <- 0.1
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

##### USER INTERFACE #####
ui <- fluidPage(
    theme=shinytheme('journal'),

    # Application title
    titlePanel("CDK2 Protein: Consensus/ML Scorings"),
    
    #fluidRow(wellPanel('Esta es una barra superior')),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            # Feature Selection and Score Type Parameters
            fluidRow(
                column(12, 
                    selectInput(
                    inputId = 'dk_score',
                    label = 'Docking score type:',
                    choices = list('Docking Score' = 'dksc',
                                   'Ligand Efficiency Score' = 'dklef') 
                ))),
            
            hr(),
            # Feature Selection and Score Type Parameters
            fluidRow(
                column(12, 
                       checkboxGroupInput(
                           inputId = 'database',
                           label = 'Ligand Database (line plot):',
                           choices = list('CSAR' = 'csar',
                                          'DUD' = 'dud',
                                          'DEKOIS2.0' = 'dekois'),
                           selected = c('csar', 'dud', 'dekois'),
                           inline = TRUE
                       ))),
            hr(),
            # Consensus/
            fluidRow(
                column(12, 
                       selectInput(
                           inputId = 'consensus',
                           label = 'Evaluation Method (CS/ML):',
                           choices = list('Consensus Scoring' = 'cons_scoring',
                                          'Manchine Learning Estimator' = 'ml_estimator'),
                       ), class = "col-md-6"),
                column(12,
                    conditionalPanel(
                        condition = "input.consensus == 'cons_scoring'",
                        checkboxGroupInput(
                           inputId = 'cs_methods',
                           label = 'CS or ML estimator used:',
                           choices = list('Rank by number' = 'rbn', 
                                          'Rank by score' = 'rbs',
                                          'Rank by rank' = 'rbr',
                                          'Rank by vote' = 'rbv_*2',
                                          'Rank by exp' = 'rexp*2'),
                           selected = c('rbn', 'rbs', 'rbr', 'rbv_*2', 'rexp*2')
                       )), 
                    conditionalPanel(
                        condition = "input.consensus == 'ml_estimator'",
                        checkboxGroupInput(
                            inputId = 'cs_methods',
                            label = 'CS or ML estimator used:',
                            choices = list('Linear SVC' = 'LinearSVC', 
                                           'RBF SVC' = 'rbfSVC',
                                           'Log. Regression' = 'LogReg',
                                           'Decision Tree' = 'Tree'),
                            selected = c('LinearSVC', 'rbfSVC', 
                                         'LogReg', 'Tree')
                        )),
                    class = "col-md-6"),
            ),
            hr(),
            
            fluidRow(
                column(12,
                  selectInput(
                    inputId = 'sel_feat_methos',
                    label = 'Method used for Conformational Selection',
                    choices = list('K-means Medoids' = 'kmeans', 
                                   'Recursive Feature Selection' = 'rfe',
                                   'Correlated Features' = 'correlated',
                                   'Random Selection' = 'random'),
                    selected = 'kmeans'
                ), class = "col-md-12")),
            fluidRow(
                column(12, radioButtons(
                    inputId = 'mds_subspace',
                    label = 'cMDS subspace (only affects k-means selection):',
                    choices = list('Pisani Resiudes' = 'pisani', 
                                   'Pocket Residues' = 'pocket'),
                    selected = 'pisani',
                    inline = TRUE
                ), class = "col-md-12")),
            width = 3
        ),

        # PLOTTING PANEL
        mainPanel(
            # Main Plot: Line plot for Consensus and ML AUC values
            fluidRow(
                div(h4('Title of the LinePlot',
                   class = 'text-center'),
                   style = 'margin-bottom: -20px; z-index: 100'),
                plotlyOutput(
                outputId = 'linearPlot')
            ),
            hr(),
            # MDS and Violin Plots to view conformations selected
            fluidRow(
                # MDS Plot
                column(12,
                    fluidRow(
                       div(h4(textOutput("mds_title"), 
                          class = 'text-center'),
                          style = 'margin-bottom: -20px; z-index: 100'),
                       plotlyOutput(outputId = 'mdsPlot')),
                       class = "col-md-6"),
                # Violin Plots of AUC Values
                column(12, 
                    fluidRow(
                       div(h4(textOutput('viol_title'),
                         class = 'text-center'),
                       style = 'margin-bottom: -20px; z-index: 100'),
                       plotlyOutput(outputId = 'swarmPlot')),
                       class = "col-md-6")
            ),
           # verbatimTextOutput("click")
           width = 9
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
    #*** Linear Plot ***
     scores_data <- reactive({
        score_ <- input$dk_score
        cs_methods <- input$cs_methods
        # Get the value from the radio button
        rd_value = input$sel_feat_methos
        # if kmeans is inside the values, update the value with mds_sub
        if('kmeans' == rd_value) {
            # Get subspace
            mds_sub <- input$mds_subspace
            # Update the value
            rd_value <- paste0('kmeans', '-', mds_sub)
        }
        # Filter the requested values
        data <- df_consensus %>%
            filter(X0 %in% c(score_) & X1 %in% c(rd_value) &
                   X2 %in% input$database & X3 %in% c(cs_methods))
        # Get the AUC values, transpose and name the columns
        data_df <- as.data.frame(t(data[,-1:-5]))
        colnames(data_df) <- apply(data, 1, get_col_names, fig = fig)
        return(data_df)
    })
     
    #*** MDS Parameters ***
    mds_subspaces <- reactive({
        paste0(input$mds_subspace, c('.dim1', '.dim2'))
    })
    
    dk_score <- reactive({
        score_ <- input$dk_score
        dk_columns_ <- grep(score_, colnames(mds_df), ignore.case =T, value = T)
    })
    
    #####################
    #***** OUTPUTS *****#
    
    #***** Text updaters
    # MDS Updater
    observe({
        mds_sub <- input$mds_subspace
        text_ <- switch(mds_sub,
                        'pisani' = 'Pisani Residues',
                        'pocket' = 'Pocket Residues')
        mds_title <- paste0('cMDS subspace (', text_, ')')
        output$mds_title <- renderText({mds_title})
    })
    # Violin plor Updater
    observe({
        dkscore <- input$dk_score
        text_ <- switch(dkscore,
                        'dksc' = 'Docking Score',
                        'dklef' = 'Ligand Efficency Score')
        viol_title <- paste0('AUC Values (Vinardo ', text_, ')')
        output$viol_title <- renderText({viol_title})
    })
    
    #***** Linear Plots
    output$linearPlot <- renderPlotly({
        
        fig <- plot_ly(type = 'scatter', mode = 'lines')
        # Draw the reference values
        ref_values_ <- ref_values_dk_scores[[input$dk_score]] # To reactive
        database_names <- names(ref_values_) # To reactive
        for(db in database_names){
            fig <- fig %>% 
                add_segments(x = 0, xend = 402, 
                             y = ref_values_[[db]], yend = ref_values_[[db]], 
                             showlegend = FALSE, name = db,  opacity = 0.5,
                             line = list(color = 'black', dash = 'dot', 
                                          linewidth = 3)) %>%
                add_text(x = 50, y = ref_values_[[db]], 
                         text = paste('max dk', toupper(db)),
                         showlegend = FALSE, textposition = "top")
        }
        # Draw the AUC results given the method
        data_df <- scores_data()
        for(col in colnames(data_df)){
            line_ <- switch (strsplit(col, ' ')[[1]][1],
                             'CSAR' = 'solid',
                             'DUD' = 'solid',
                             'DEKOIS' = 'solid')
            fig <- fig %>% add_trace(y = data_df[[col]], 
                            name = col, line = list(dash = line_))
        }
        fig %>% layout(xaxis = ax_lp, yaxis = yax_lp, 
                       paper_bgcolor='rgba(0,0,0,0)', 
                       legend = list(title = 
                                list(text = '<b>Database & Method:</b>'))) %>%
            config(modeBarButtonsToRemove = conf_, displaylogo = FALSE) 
    })
    
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
                   paper_bgcolor='rgba(0,0,0,0)', 
                   legend = list(title = 
                            list(text = '<b>Protein<br>Conformations:</b>'))) %>%
                   # title = list(text = paste0('cMDS subspace (',
                   #                            input$mds_subspace, ')'), x = 0.1)) %>%
            config(modeBarButtonsToRemove = conf_, displaylogo = FALSE) %>%
            event_register("plotly_click")
        
    })
    
    #*** Observe Events
    observeEvent(input$consensus, {
        if(input$consensus == 'cons_scoring'){
            # Show Consensus Metrics
            shinyjs::show(id = 'cs_methods')
        } else {
            # Show Machine Learning Estiamtors
            shinyjs::hide(id = 'cs_methods')
        }
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
                             side = 'positive',
                             text = text)
        i = 0
        for (score_ in score_cols) {
            i = i + 1
            fig_swarm <- add_trace(fig_swarm, y = mds_df[[score_]], pointpos = -1, 
                      marker = list(size = marker_size), 
                      name = strsplit(score_, '_')[[1]][1],
                      color = I(pal_violin[i]),
                      box = list(visible = T), key = rownames(mds_df))
        }
        fig_swarm <- layout(fig_swarm, xaxis = xax_sw, yaxis = yax_sw,
                            paper_bgcolor='rgba(0,0,0,0)', 
                            dragmode = 'select', 
                            legend = list(title = 
                                    list(text = '<b>Ligand<br>Databases:</b>'))) %>% 
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
                                   hoverinfo='skip',
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
