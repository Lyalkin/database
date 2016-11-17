###
# This file contains the functions to compute the potential distribution of species, according to climatic variables
# The functions provided are:
#   - copmute the potential area according to a combination of climatic variables (/potential_with/ and /potential_without/)
#   - Compute the Shapley inclusion score for one variable (/shapley_score_with/ and /shapley_score_without/)
# These functions can handle the model without and with interactions, according to their suffixes:
#   - the suffix "_with" appended to the function name, for the model with interactions
#   - "_without", for the model without interaction
# Also, a function is provided to output a color bar plot (/barplot_shapley/), possibly with computing bootstrapped confidence intervals
###


# This computes the potential area based on the model WITH interactions
# Args: vars_selected, a vector of length 1:19 representing the climatic variables used
#       fia_data, a reference class with the data (required for quicker processing)
# Output: the prediction vector (binary vector with as many items as there are plots in the data)
potential_with = function (vars_selected, clim_data)
{
  n_vars_selected = length(vars_selected)
  dimensions = paste0('Bio', vars_selected, 's_bin')
  names(dimensions) = paste0('dim', 1:n_vars_selected)
  fia_db_with_Bins_marker_T = rep(0, clim_data$n_T)
  fia_db_with_Bins_marker_all = rep(0, clim_data$n_all)
  for (i in 1:n_vars_selected) {
    if (i <= 10) {
      fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + clim_data$fia_db_with_Bins_T[,dimensions[i]]*res^(i-1)
      fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + clim_data$fia_db_with_Bins[,dimensions[i]]*res^(i-1) 
    } else {
      fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + clim_data$fia_db_with_Bins_T[,dimensions[i]]/res^(i-10)
      fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + clim_data$fia_db_with_Bins[,dimensions[i]]/res^(i-10) 
    }
  }
  prediction = fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T))
  return(prediction)
}



# This computes the potential area based on the model WITHOUT interactions
# Args: vars_selected, a vector of length 1:19 representing the climatic variables used
#       fia_data, a reference class with the data (required for quicker processing)
# Output: the prediction vector (binary vector with as many items as there are plots in the data)
potential_without = function (vars_selected, clim_data)
{
  n_vars_selected = length(vars_selected)
  dimensions = paste0('Bio', vars_selected, 's_bin')
  names(dimensions) = paste0('dim', 1:n_vars_selected)
  fia_db_with_Bins_marker_all = rep(0, clim_data$n_all)
  for (i in 1:n_vars_selected) {
    potential_climates_i = unique(clim_data$fia_db_with_Bins_T[,dimensions[i]])
    potential_dist = which(clim_data$fia_db_with_Bins[,dimensions[i]] %in% potential_climates_i)
    fia_db_with_Bins_marker_all[potential_dist] = fia_db_with_Bins_marker_all[potential_dist] + 1
  }
  prediction = (fia_db_with_Bins_marker_all == n_vars_selected)
  return(prediction)
}




# This implements the Shapley inclusion score based on the model WITH interactions
# Args: vars_selected, a vector of length 1:18 - they are indices of climatic variables
#       var_toggle, an extra variable for which we compute the Shapley inclusion score
#       fia_data, a reference class with the data (required for quicker processing)
# Output: the shapley inclusion score of 'var_toggle' to 'vars_selected', according to 'data'
shapley_score_with = function (vars_selected, var_toggle, clim_data)
{
  n_vars_selected = length(vars_selected)
  all_vars = c(vars_selected, var)
  dimensions = paste0('Bio', all_vars, 's_bin')
  names(dimensions) = paste0('dim', 1:(n_vars_selected+1))
  fia_db_with_Bins_marker_T = rep(0, clim_data$n_T)
  fia_db_with_Bins_marker_all = rep(0, clim_data$n_all)
  # compute the score without the additional var:
  for (i in 1:n_vars_selected) {
    if (i <= 10) {
      fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + clim_data$fia_db_with_Bins_T[,dimensions[i]]*res^(i-1)
      fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + clim_data$fia_db_with_Bins[,dimensions[i]]*res^(i-1) 
    } else {
      fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + clim_data$fia_db_with_Bins_T[,dimensions[i]]/res^(i-10)
      fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + clim_data$fia_db_with_Bins[,dimensions[i]]/res^(i-10) 
    }
  }
  score_without = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
  # compute the score with the additional var:
  i = i + 1
  if (i <= 10) {
    fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + clim_data$fia_db_with_Bins_T[,dimensions[i]]*res^(i-1)
    fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + clim_data$fia_db_with_Bins[,dimensions[i]]*res^(i-1) 
  } else {
    fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + clim_data$fia_db_with_Bins_T[,dimensions[i]]/res^(i-10)
    fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + clim_data$fia_db_with_Bins[,dimensions[i]]/res^(i-10)
  }
  score_with = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
  # shapley value for the difference:
  shapleys = (score_without - score_with) * 
    ( factorial(n_vars_selected) * factorial((clim_data$variables_nb-1)-n_vars_selected) / factorial(clim_data$variables_nb) ) 
  return(shapleys)
}



# This implements the Shapley inclusion score based on the model WITHOUT interactions
# Args: vars_selected, a vector of length 1:18 - they are indices of climatic variables
#       var_toggle, an extra variable for which we compute the Shapley inclusion score
#       fia_data, a reference class with the data (required for quicker processing)
# Output: the shapley inclusion score of 'var_toggle' to 'vars_selected', according to 'data'
shapley_score_without = function (vars_selected, var_toggle, clim_data)
{
  n_vars_selected = length(vars_selected)
  all_vars = c(vars_selected, var)
  dimensions = paste0('Bio', all_vars, 's_bin')
  names(dimensions) = paste0('dim', 1:(n_vars_selected+1))
  fia_db_with_Bins_marker_all = rep(0, clim_data$n_all)
  # compute the score without the additional var:
  for (i in 1:n_vars_selected) {
    potential_climates_i = unique(clim_data$fia_db_with_Bins_T[,dimensions[i]])
    potential_dist = which(clim_data$fia_db_with_Bins[,dimensions[i]] %in% potential_climates_i)
    fia_db_with_Bins_marker_all[potential_dist] = fia_db_with_Bins_marker_all[potential_dist] + 1
  }
  score_without = sum(fia_db_with_Bins_marker_all == n_vars_selected)
  # compute the score with the additional var:
  i = i + 1
  potential_climates_var = unique(clim_data$fia_db_with_Bins_T[,dimensions[i]])
  potential_dist = which(clim_data$fia_db_with_Bins[,dimensions[i]] %in% potential_climates_i)
  fia_db_with_Bins_marker_all[potential_dist] = fia_db_with_Bins_marker_all[potential_dist] + 1
  score_with = sum(fia_db_with_Bins_marker_all == (n_vars_selected+1))
  # shapley value for the difference:
  shapleys = (score_without - score_with) * 
    ( factorial(n_vars_selected) * factorial((clim_data$variables_nb-1)-n_vars_selected) / factorial(clim_data$variables_nb) ) 
  return(shapleys)
}



# This output a barplot for the Shapley values.
# Args: shapleys, a matrix of shapley inclusion scores
#       bootstrap_n, the optional number of bootstrap iterations for the confidence intervals (default: 0, i.e. no confidence interval is computed) 
#       bootstrap_alpha, the confidence level (default: 95%)
# Output: the shapley inclusion score of 'var_toggle' to 'vars_selected', according to 'data'
barplot_shapley = function(shapleys, bootstrap_n = 0, bootstrap_conf = 0.95)
{
  require(boot) ## If error, install: install.packages("boot")
  
  n_done = which(is.na(shapleys[,ncol(shapleys)]))[1] - 1
  if (is.na(n_done)) {
    n_done = nrow(shapleys)
  }
  
  # this variable encodes whether the bioclimatic variable is linked to temperature (1), precipitation (2) or if it is some measure of variability (3)
  bioclim_types = c(1, 3, 3, 3, 1, 1, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2)[1:ncol(shapleys)]
  # the color of the different types:
  colormap = c('orange','lightblue', 'purple')
  # the variable names written in a friendly way
  all_variable_names = c('Annual Mean Temperature',
                         'Mean Diurnal Range',
                         'Isothermality',
                         'Temperature Seasonality',
                         'Max Temperature of Warmest Month',
                         'Min Temperature of Coldest Month',
                         'Temperature Annual Range',
                         'Mean Temperature of Wettest Quarter',
                         'Mean Temperature of Driest Quarter',
                         'Mean Temperature of Warmest Quarter',
                         'Mean Temperature of Coldest Quarter',
                         'Annual Precipitation',
                         'Precipitation of Wettest Month',
                         'Precipitation of Driest Month',
                         'Precipitation Seasonality',
                         'Precipitation of Wettest Quarter',
                         'Precipitation of Driest Quarter',
                         'Precipitation of Warmest Quarter',
                         'Precipitation of Coldest Quarter')[1:ncol(shapleys)]
  
  if (bootstrap_n == 0) {
    
    # simple version, no bootstrap-derived confidence intervals
    scores_exact = colMeans(shapleys, na.rm=T)
    order_exact = order(scores_exact, decreasing = T)
    bp = barplot(scores_exact[order_exact], col=colormap[bioclim_types[order_exact]], names.arg=NA, las=1, ylab='Shapley value')
    text(bp, 0, paste0(' ',all_variable_names[order_exact]), srt=90, adj=0, xpd=T, font=3)
    
  } else {
    
    # compute the confidence interval
    get_mean_subsample = function(data, indices) {
      return(colMeans(data[indices, ]))
    }
    results = boot(data=shapleys[1:n_done, ], statistic=get_mean_subsample, R=bootstrap_n)
    confs = sapply(1:ncol(shapleys), function(i){
      tmp=boot.ci(results, type="norm", index=i, conf=bootstrap_conf)
      return(c(tmp$t0, tmp$normal))}
    )
    confs = rbind(confs, 1:ncol(shapleys))
    row.names(confs) = c('statistic', 'conf_int', 'min_bound', 'max_bound', 'climatic_var_i')
    colnames(confs) = all_variable_names
    ordered_shap = confs[, order(confs[1,], decreasing = T)]
    scores = unlist(ordered_shap[1,])
    
    bp = barplot(scores, col=colormap[bioclim_types[unlist(ordered_shap[5,])]], names.arg=NA, las=1, ylab='Shapley value')
    text(bp, 0, paste0(' ',names(ordered_shap[1,])), srt=90, adj=0, xpd=T, font=3)
    require('Hmisc')
    errbar(bp, unlist(ordered_shap[1,]), unlist(ordered_shap[3,]), unlist(ordered_shap[4,]), add=T, xpd=T, pch=NA)
  }
  
  legend('topright', fill=colormap, c('Temperature-related', 'Precipitation-related', 'Variability-related'), bty='n', xpd=T)
}

