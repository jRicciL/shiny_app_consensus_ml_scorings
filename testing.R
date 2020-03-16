# plot the cMDS subspace with plotly
library(ggplot2)
library(plotly)
library(dplyr)

##### CONSENSUS SCORIGS
df_consensus <- read.csv('data/kmeans_cosensus_scorings.csv')



##### MDS
# Load the conformatiosn data
# Pdb_id, label conf, pk volume, 
df_prot <- read.csv('data/df_pdb_cdk2_features_and_scores.csv', header = T,
                    na.strings=c("","NA"), stringsAsFactors =  F)
# replacing NA values
df_prot[is.na(df_prot)] <- 0
# Numeric columns
df_prot$Inhib_mass <- as.numeric(df_prot$Inhib_mass)^3

# Load the mds data
mds_data <- data.frame(t(read.csv('data/mds_subspaces_pis_pk.csv', header = FALSE)))
subspace_names = c('pisani', 'pocket')
n_dims = paste0('dim', c(1:2))
# Name the columns dataframe
colnames(mds_data) <- paste(rep(subspace_names, each = length(n_dims)), n_dims, sep = '.')

# Combine in a unique dataframe
df_plot <- cbind(mds_data, df_prot)
rownames(df_prot) <- df_prot$X



# Need to save this dataframe as an R object
list_objs = list(mds_df = df_plot)
saveRDS(list_objs, 'data/data.rds')

# colors
pal <- c("#f64b3c", "#5b8c5a", "#0f4c75", "#feb72b", "#29c7ac")
pal <- setNames(pal, c('active', 'inact_a', 'inact_b', 'inact_ope', 'dfg_out'))



ax_mds <- list(linecolor = toRGB('black'), automargin = T, 
               scaleratio =1, scaleanchor = "x",  dtick = 0.25, 
               zerolinecolor = toRGB('grey'), zerolinewidth = 2, mirror = T, showline = T)
yax_mds <- ax_mds; yax_mds[['title']] <- 'Second Dimension'
xax_mds <- ax_mds; xax_mds[['title']] <- 'First Dimension'

yax_sw <- ax_mds; yax_sw[['title']] <- 'ROC-AUC'; yax_sw[['range']] <- c(0.3, 1)
xax_sw <- ax_mds; xax_sw[['title']] <- 'Databases'

points <- c(1:4)
sliced_df <- df_plot[points, ]

fig_mds <- plot_ly(type = 'scatter', mode = 'markers', 
                   data = df_plot, colors = pal) %>% 
  add_trace(x = df_plot[['pocket.dim1']], y = df_plot[['pocket.dim2']],
               color = df_plot[['Labels_conf']], 
               text = ~paste('<b>PDB id:</b>', df_plot$X,
                             '<br><b>Conf:</b>', df_plot$Labels_conf,
                             '<br><b>Ligand:</b>', df_plot$Inhib),
               size = 0.1,
               hovertemplate = paste('%{text}')) %>% 
  # add_trace(x = sliced_df[['pocket.dim1']], y = sliced_df[['pocket.dim2']],
  #           marker = list(color = 'rgba(0, 0, 0, 0)',
  #                         line = list(color = 'rgba(0, 0, 0, 1)',
  #                                     width = 2))) %>%
  layout(xaxis = xax_mds, yaxis = yax_mds, dragmode = 'pan')
fig_mds



### BOX PLOT
text = ~paste('<b>PDB id:</b>', df_plot$X,
              '<br><b>Conf:</b>', df_plot$Labels_conf,
              '<br><b>Ligand:</b>', df_plot$Inhib)




marker_size = 5
fig_swarm <- plot_ly(type = 'violin', data = df_plot, 
                     points = "all", jitter = 1, colors = pal,
                     selectedpoints = c(1,2), text = text,
                     box = list(visible = T)) %>%
  add_trace(y = df_plot[['CSAR_vrd_dksc']], pointpos = 0, 
            marker = list(size = marker_size), fillcolor = I("#8dd3c7")) %>%
  add_trace(y = df_plot[['DUD_vrd_dksc']], pointpos = 0, 
            marker = list(size = marker_size)) %>%
  add_trace(y = df_plot[['DEKOIS_vrd_dksc']], pointpos = 0, 
            marker = list(size = marker_size)) %>% 
  layout(xaxis = xax_sw, yaxis = yax_sw, dragmode = 'pan', showlegend = FALSE) %>% 
  highlight(on = "plotly_hover", off = "plotly_doubleclick") %>%
  config(modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "select2d", "toImage", "autoScale2d",
                                    "toggleSpikelines", "hoverCompareCartesian", "hoverClosestCartesian"),
         displaylogo = FALSE)
fig_swarm


subplot(fig_mds, fig_swarm, nrows = 2)  %>%
  highlight(
    dynamic = TRUE, 
    selectize = TRUE, 
    selected = attrs_selected(opacity = 0.3)
  )


df_plot %>%
  add_markers(x = ~jitter(as.numeric(group)), y = ~xval, color = ~group,
              marker = list(size = 6),
              hoverinfo = "text",
              text = ~paste0("Group: ",group,
                             "<br>xval: ",xval),
              showlegend = FALSE) 


dat <- data.frame(xval = sample(100,1000,replace = TRUE),
                  group = as.factor(sample(c("a","b","c"),1000,replace = TRUE)))

dat %>%
  plot_ly() %>% 
  add_markers(x = ~jitter(as.numeric(group)), y = ~xval, color = ~group,
              marker = list(size = 6),
              hoverinfo = "text",
              text = ~paste0("Group: ",group,
                             "<br>xval: ",xval),
              showlegend = FALSE) %>% 
  layout(legend = list(orientation = "h",
                       x =0.5, xanchor = "center",
                       y = 1, yanchor = "bottom"
  ),
  xaxis = list(title = "Group",
               showticklabels = FALSE))


###
data(txhousing, package = "ggplot2")

# declare `city` as the SQL 'query by' column
tx <- highlight_key(txhousing, ~city)

# initiate a plotly object
base <- plot_ly(tx, color = I("black")) %>% 
  group_by(city)

hist <- add_histogram(
  base,
  x = ~median, 
  histnorm = "probability density"
)
subplot(time_series, hist, nrows = 2) %>%
  layout(barmode = "overlay", showlegend = FALSE) %>%
  highlight(
    dynamic = TRUE, 
    selectize = TRUE, 
    selected = attrs_selected(opacity = 0.3)
  )
