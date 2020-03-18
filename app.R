library(shiny)
library(shinythemes)
library(plotly)
library(dplyr)
library(shinyjs)
library(shinycssloaders)

options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=1)

# Load the date once
# Load the RDS file
list_objs <- readRDS('data/data.rds')

# Get the consensus and ml dataframes
df_consensus <- list_objs$df_consensus
list_ml <- list_objs$list_ml

selected_features <- list_objs$selected_features

# List of database names
db_names <- list('CSAR' = 'csar', 'DUD' = 'dud', 'DEKOIS2.0' = 'dekois')

# Reference values
ref_values_dk_scores <- list(
    'dksc'  = list('csar' = 0.8486, 'dud' = 0.675, 'dekois' = 0.7904),
    'dklef' = list('csar' = 0.9287, 'dud' = 0.7788, 'dekois' = 0.7342))

# funtion to use apply
get_col_names <- function(x, fig){
    trace_name <- paste0(c(toupper(x[5]), x[6]), collapse = ' ')
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
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),    

    # Application title
    titlePanel("CDK2 Protein: Consensus/ML Scorings"),
    
    #fluidRow(wellPanel('Esta es una barra superior')),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            # Feature Selection and Score Type Parameters
            h3('Docking score type'),
            fluidRow(
                column(12, 
                    selectInput(
                    inputId = 'dk_score',
                    label = 'Select the Docking Score type:',
                    choices = list('Docking Score' = 'dksc',
                                   'Ligand Efficiency Score' = 'dklef') 
                ))),
            hr(),
            h3('Evaluation method'),
            fluidRow(
                column(12, 
                       selectInput(
                           inputId = 'consensus',
                           label = 'Evaluation Method (CS/ML):',
                           choices = list('Manchine Learning' = 'ml_estimator',
                                          'Consensus Scoring' = 'cons_scoring'
                                          ),
                       ), class = "col-md-6"),
                column(12,
                    conditionalPanel(
                        condition = "input.consensus == 'cons_scoring'",
                        checkboxGroupInput(
                           inputId = 'cs_methods',
                           label = 'CS or ML estimator used:',
                           choices = list('Rank by number' = 'rbn', 
                                          'Rank by score'  = 'rbs',
                                          'Rank by rank'   = 'rbr',
                                          'Rank by vote'   = 'rbv_*2',
                                          'Rank by exp'    = 'rexp*2'),
                           selected = c('rbn', 'rbs', 'rbr', 'rbv_*2', 'rexp*2')
                       )), 
                    conditionalPanel(
                        condition = "input.consensus == 'ml_estimator'",
                        checkboxGroupInput(
                            inputId = 'ml_methods',
                            label = 'CS or ML estimator used:',
                            choices = list('Linear SVC' = 'LinearSVC', 
                                           'RBF SVC' = 'rbfSVC',
                                           'Log. Regression' = 'LogRg'),
                            selected = c('LinearSVC', 'rbfSVC', 
                                         'LogRg')
                        )),
                    class = "col-md-6"),
            ),
            hr(),
            # Feature Selection and Score Type Parameters
            h3('Ligand databases'),
            fluidRow(
                column(12, 
                    conditionalPanel(
                       condition = "input.consensus == 'cons_scoring'",
                       checkboxGroupInput(
                           inputId = 'database',
                           label = 'Ligand Database (Consensus Scoring):',
                           choices = list('CSAR' = 'csar',
                                          'DUD' = 'dud',
                                          'DEKOIS2.0' = 'dekois'),
                           selected = c('csar', 'dud', 'dekois'),
                           inline = TRUE
                       )
                    ),
                    conditionalPanel(
                      condition = "input.consensus == 'ml_estimator'",
                      fluidRow(
                        column(12,
                          selectInput(
                            inputId = 'train_db',
                            label = 'Training Database:',
                            choices = list('CSAR' = 'csar',
                                           'DUD' = 'dud',
                                           'DEKOIS2.0' = 'dekois'),
                          ),
                          class = "col-md-6"
                        ),
                        column(12,
                          checkboxGroupInput(
                            inputId = 'test_db',
                            label = 'Testing Databases:',
                            choices = db_names
                          ),
                          class = "col-md-6"
                        )
                      )
                    )
                  )
                ),
            hr(),
            h3('Conformational Selection'),
            fluidRow(
                column(12,
                  selectInput(
                    inputId = 'sel_feat_methods',
                    label = 'Method used for Conformational Selection',
                    choices = list('K-means Medoids' = 'kmeans',
                                   'Recursive Feature Selection' = 'rfe',
                                   'Correlated Features' = 'correlated',
                                   'Random Selection (means)' = 'random'),
                    selected = 'kmeans'
                ), class = "col-md-12")),
            conditionalPanel(
              condition = "input.sel_feat_methods == 'kmeans'",
              fluidRow(
                  column(12, radioButtons(
                      inputId = 'mds_subspace',
                      label = 'cMDS subspace (only affects k-means selection):',
                      choices = list('Pisani Resiudes' = 'pisani', 
                                     'Pocket Residues' = 'pocket'),
                      selected = 'pisani',
                      inline = TRUE
                  ), class = "col-md-12")),
            ),
            fluidRow(
              column(12,
                     sliderInput(
                       inputId = 'feature_slider',
                       label = 'Number of features used for the evaluation:',
                       min = 0, max = 200,
                       value = 1
                     ), class = "col-md-12")),
            
            width = 3
        ),

        # PLOTTING PANEL
        mainPanel(
            # Main Plot: Line plot for Consensus and ML AUC values
            fluidRow(
                div(h4(textOutput("linePlot_title"),
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
                       withSpinner(plotlyOutput(outputId = 'mdsPlot'))),
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
    ),
    div(class = "footer",
        'J. Ricci-Lopez, 2020'
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
  
    #*** Linear Plot ***
     scores_data <- reactive({
        score_ <- input$dk_score
        
        # Get the value from the radio button
        rd_value = input$sel_feat_methods
        # if kmeans is inside the values, update the value with mds_sub
        if('kmeans' == rd_value) {
            # Get subspace
            mds_sub <- input$mds_subspace
            # Update the value
            rd_value <- paste0('kmeans', '-', mds_sub)
        }
        # Select the correct methods
        methods <- switch (input$consensus,
            'cons_scoring' =  input$cs_methods,
            'ml_estimator' = input$ml_methods
        )
        
        # Selects the database
        databases <- switch (input$consensus,
             'cons_scoring' = input$database,
             'ml_estimator' = input$test_db
        )
         
        # Select the dataset
        main_df <- switch (input$consensus,
            'cons_scoring' = df_consensus,
            'ml_estimator' = list_ml[[input$train_db]]
        )
        
        # Filter the requested values
        data <- main_df %>%
            filter(X1 %in% c(score_) & X2 %in% c(rd_value) &
                   X3 %in% c(databases) & X4 %in% c(methods))
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
    # Line plot tile
    observe({
      # Docking Score
      dk_score <- switch(input$dk_score,
          'dksc' = '(DkSc)',
          'dklef' = '(DkLEf)')
      # Evaluation
      eval_method <- switch (input$consensus,
          'ml_estimator' = 'Machine Learning',
          'cons_scoring' = 'Consensus Scoring'
      )
      # Training if ML
      train_db <- ''
      if(input$consensus == 'ml_estimator'){
        train_db <- paste0('[', toupper(input$train_db), ' training]')
      }
      # Feature Selection
      feature_sel <- switch (input$sel_feat_methods,
           'kmeans' = 'K-means Selection',
           'rfe' = 'Recursive Feature Selection',
           'correlated' = 'Dropping Correlated Features',
           'random' = 'Random Selection'
      )
      linePlot_title <- paste(eval_method, train_db, dk_score, '-', feature_sel)
      output$linePlot_title <- renderText({linePlot_title})
    })
    
    #***** Linear Plots
    output$linearPlot <- renderPlotly({
        
        fig <- plot_ly(type = 'scatter', mode = 'lines')
        # Draw the reference values
        ref_values_ <- ref_values_dk_scores[[input$dk_score]] # To reactive
        #database_names <- names(ref_values_) # To reactive
        # Selects the databases
        database_names <- switch (input$consensus,
             'cons_scoring' = input$database,
             'ml_estimator' = input$test_db
        )
        
        for(db in database_names){
            fig <- fig %>% 
                add_segments(x = 0, xend = 402, 
                             y = ref_values_[[db]], yend = ref_values_[[db]], 
                             showlegend = FALSE, name = db,  opacity = 0.5,
                             line = list(color = 'black', dash = 'dash', 
                                          linewidth = 3)) %>%
                add_text(x = 50, y = ref_values_[[db]], 
                         text = paste('max dk', toupper(db)),
                         showlegend = FALSE, textposition = "top")
        }
        # Draw the AUC results given the method
        data_df <- scores_data()
        for(col in colnames(data_df)){
            line_ <- switch (strsplit(col, ' ')[[1]][1],
                             'CSAR'   = 'solid',
                             'DUD'    = 'solid',
                             'DEKOIS' = 'solid')
            fig <- fig %>% add_trace(y = data_df[[col]], 
                            name = col, line = list(dash = line_)
                           )
        }
        x_pos <- input$feature_slider
        fig <- fig %>% 
          add_segments(x = x_pos, xend = x_pos,
                       y = 0, yend = 1.0, opacity = 0.8,
                       showlegend = FALSE,
                       line = list(color = '#f64a3a', 
                                   dash = 'dot', linewidth = 6))
        fig %>% layout(xaxis = ax_lp, yaxis = yax_lp, 
                       paper_bgcolor = 'rgba(0,0,0,0)', 
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
            add_trace(x = NA, y = NA, name = 'Selected') %>% # Empty trace
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
    
    # Updates the selected
    
    # Oserve Event to uodated Databases for training
    observeEvent(input$train_db, {
      db_choices <- db_names
      db_choices[db_choices == input$train_db] <- NULL
      updateCheckboxGroupInput(session = session,
                               inputId = 'test_db',
                               choices = db_choices,
                               selected = db_choices)
    })
    
    
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
        points <- selected_poins()
        fig_swarm <- plot_ly(type = 'violin', data = mds_df, 
                             points = "all", jitter = 1,
                             side = 'positive',
                             text = text, selectedpoints = points)
        i = 0
        for (score_ in score_cols) {
            i = i + 1
            fig_swarm <- add_trace(fig_swarm, y = mds_df[[score_]], pointpos = -1, 
                      marker = list(size = marker_size), 
                      name = strsplit(score_, '_')[[1]][1],
                      color = I(pal_violin[i]),
                      text = text_hover,
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
    
    # Reactive selection of points
    selected_poins <- reactive({
      n_feat <- input$feature_slider
      
      method <- input$sel_feat_methods # Get the method
      if('kmeans' == method) {
        # Get subspace
        mds_sub <- input$mds_subspace
        # Update the value
        method <- paste0('kmeans', '-', mds_sub)
      }
      features <- selected_features[[method]][[n_feat]] + 1
      features
    })
    
    
    # Selection from slider
    observeEvent({input$feature_slider
                  input$sel_feat_methods
                  input$mds_subspace
                  }, {
      # Get the correct features
      n_feat <- input$feature_slider
      
      method <- input$sel_feat_methods # Get the method
      if('kmeans' == method) {
        # Get subspace
        mds_sub <- input$mds_subspace
        # Update the value
        method <- paste0('kmeans', '-', mds_sub)
      }
      
      mds_ <- mds_subspaces()
      features <- selected_features[[method]][[n_feat]] + 1
      sliced_df <- mds_df[features, ]
      # # Remove the previous selected point if it exist
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
