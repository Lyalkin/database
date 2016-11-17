# set the working directory
setwd("~/tmptodelete/model_without_interactions/") # change this on your computer
  
# source the functions
source('model_fcts.R')

# Create a reference class to encapsulate the data (enables passing large objects by reference)
clim_data_generator = setRefClass("clim_data_class", fields = list(
  fia_db_with_Bins = "data.frame", 
  fia_db_with_Bins_T = "data.frame",
  n_all = "numeric", 
  n_T = "numeric", 
  res = "numeric",
  variables_nb = "numeric"))



########################################################################################
# TEST OF THE POTENTIAL DISTRIBUTION PREDICTION OF THE MODEL WITH/WITHOUT INTERACTIONS #
########################################################################################

# Test with an example that looks like this, with 4 data points (a,b,c,d) coded along DIM 1 and DIM 2
#           DIM 1
#          1 2 3 4
#       #############
# D  1  #  . . . d  #
# I  2  #  . . . .  #
# M  3  #  . b . .  #
#    4  #  a c . .  #
# 2     #############
#

TEST = expand.grid(Bio1s_bin=1:4, Bio2s_bin=1:4)
TEST$EXAMPLE = F
TEST$EXAMPLE[c(13, 10, 14, 4)] = T

# Generate the fake climatic data for this example
clim_data = clim_data_generator(
  fia_db_with_Bins = TEST,
  fia_db_with_Bins_T = TEST[which(TEST[, 'EXAMPLE']==T), ],
  n_all = nrow(TEST),
  n_T = length(which(TEST[, 'EXAMPLE']==T)),
  res = 4,
  variables_nb = 2
)

# Get the prediction on the model with interactions:
prediction_with = potential_with(c(1,2), clim_data)
# Plot it on the commandline:
cc=1
for (dim1 in 1:4) {
  for (dim2 in 1:4) {
    if (prediction_with[cc])
      cat('x ')
    else
      cat('. ')
    cc = cc+1
  }
  cat('\n')
}
### OUTPUT:
# . . . x 
# . . . . 
# . x . . 
# x x . . 


# Get the prediction on the model without interaction:
prediction_without = potential_without(c(1,2), clim_data)
# Plot it on the commandline:
cc=1
for (dim1 in 1:4) {
  for (dim2 in 1:4) {
    if (prediction_without[cc])
      cat('x ')
    else
      cat('. ')
    cc = cc+1
  }
  cat('\n')
}
### OUTPUT:
# x x . x 
# . . . . 
# x x . x 
# x x . x 



################################################################################################
# TEST OF THE MODEL WITH/WITHOUT INTERACTIONS, COMPUTED WITH THE EXHAUSTIVE/RESAMPLED APPROACH #
################################################################################################


# load the example data
load('example.RData')

# Get a subset of the climatic data for Species 98
fia_db_with_Bins = fia_db_with_Bins[sample.int(nrow(fia_db_with_Bins),50000),] # keep only a subset, for faster testing
clim_data = clim_data_generator(
  fia_db_with_Bins = fia_db_with_Bins,
  fia_db_with_Bins_T = fia_db_with_Bins[which(fia_db_with_Bins[, data_name]==T), ],
  n_all = nrow(fia_db_with_Bins),
  n_T = length(which(fia_db_with_Bins[, data_name]==T)),
  res = res,
  variables_nb = 10 # keep only a few variables, for faster testing (should be 19 for all climatic variables)
)


# Resampling approach to approximate the shapley values, model WITH interactions
shapleys_WITH_resampled = matrix(NA, ncol=clim_data$variables_nb, nrow=5000)
pb = txtProgressBar(max=nrow(shapleys_WITH_resampled), style = 3)
for (trial in 1:nrow(shapleys_WITH_resampled)) {
  setTxtProgressBar(pb, trial)
  for (var in 1:ncol(shapleys_WITH_resampled)) {
    n_vars = sample.int((clim_data$variables_nb-1), 1)
    vars_selected = sample((1:clim_data$variables_nb)[-var], n_vars)
    shapleys_WITH_resampled[trial, var] = shapley_score_with(vars_selected, var, clim_data)
  }
}


# Resampling approach to approximate the shapley values, model WITHOUT interactions
shapleys_WITHOUT_resampled = matrix(NA, ncol=clim_data$variables_nb, nrow=5000)
pb = txtProgressBar(max=nrow(shapleys_WITHOUT_resampled), style = 3)
for (trial in 1:nrow(shapleys_WITHOUT_resampled)) {
  setTxtProgressBar(pb, trial)
  for (var in 1:ncol(shapleys_WITHOUT_resampled)) {
    n_vars = sample.int((clim_data$variables_nb-1), 1)
    vars_selected = sample((1:clim_data$variables_nb)[-var], n_vars)
    shapleys_WITHOUT_resampled[trial, var] = shapley_score_without(vars_selected, var, clim_data)
  }
}


# Exhaustive approach to compute the shapley values:
possible_comb = as.matrix(expand.grid(as.data.frame(sapply(1:(clim_data$variables_nb-1), function(...)c(T,F))))) # all the T/F combinations of 18 variables (19-1 variables)
shapleys_WITH_exact = matrix(NA, ncol=clim_data$variables_nb, nrow=nrow(possible_comb)-1)
pb = txtProgressBar(max=nrow(shapleys_WITH_exact), style = 3)
for (trial in 1:nrow(shapleys_WITH_exact)) {
  setTxtProgressBar(pb, trial)
  for (var in 1:ncol(shapleys_WITH_exact)) {
    vars_selected = which(possible_comb[trial,])
    vars_selected[which(vars_selected>=var)] = vars_selected[which(vars_selected>=var)] + 1
    shapleys_WITH_exact[trial, var] = shapley_score_with(vars_selected, var, clim_data)
  }
}


# Exhaustive approach to compute the shapley values:
possible_comb = as.matrix(expand.grid(as.data.frame(sapply(1:(clim_data$variables_nb-1), function(...)c(T,F))))) # all the T/F combinations of 18 variables (19-1 variables)
shapleys_WITHOUT_exact = matrix(NA, ncol=clim_data$variables_nb, nrow=nrow(possible_comb)-1)
pb = txtProgressBar(max=nrow(shapleys_WITHOUT_exact), style = 3)
for (trial in 1:nrow(shapleys_WITHOUT_exact)) {
  setTxtProgressBar(pb, trial)
  for (var in 1:ncol(shapleys_WITHOUT_exact)) {
    vars_selected = which(possible_comb[trial,])
    vars_selected[which(vars_selected>=var)] = vars_selected[which(vars_selected>=var)] + 1
    shapleys_WITHOUT_exact[trial, var] = shapley_score_without(vars_selected, var, clim_data)
  }
}


pdf('comparison.pdf', height=12, width=12)
par(mfrow=c(2,2))
barplot_shapley(shapleys_WITH_resampled, bootstrap_n = 1000)
title('Model WITH interaction, resampled approximation')

barplot_shapley(shapleys_WITHOUT_resampled, bootstrap_n = 1000)
title('Model WITHOUT interaction, resampled approximation')

barplot_shapley(shapleys_WITH_exact)
title('Model WITH interactions, exhaustive combinations')

barplot_shapley(shapleys_WITHOUT_exact)
title('Model WITHOUT interaction, exhaustive combinations')
dev.off()