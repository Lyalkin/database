# set the working directory
setwd("~/tmptodelete/model_without_interactions/") # change this on your computer

code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)


require('maps')
require('raster')
require('sp')
require('rgdal')


# 125 is for Red pines (Pinus resinosa)
species_number = 125


# extraction of all the plots in the FIA DB and computation of the relative BA of the species
# source('fia_extract.R')
# fia_db = fia_extract(p=data_path, species_number)
# saveRDS(fia_db, file = paste0(code_path, '/fia_db_', species_number, '.rds'))
fia_db = readRDS(file = paste0(code_path, '/fia_db_', species_number, '.rds'))


# Extraction of the climatic variables from Worldclim
#source('clim_extract.R')
#tmpppt = extract_tmpppt (fia_db, p=data_path) 
#saveRDS(tmpppt, file = paste0(code_path, '/tmpppt.rds'))
tmpppt = readRDS(file = paste0(code_path, '/tmpppt.rds'))



tmpppt_df = as.data.frame(tmpppt)


# Now make tmpppt_df_binned from tmpppt_df
source('extract_breaks_bins.R')
breaks = extract_breaks(tmpppt_df)
bins = extract_bins(tmpppt_df, breaks)

#tmpppt_df_binned is tmpppt_df where each column is binned
tmpppt_df_binned = as.data.frame(bins[[1]])
names(tmpppt_df_binned) = (paste0(names(tmpppt_df)[1],'_bin',collapse = ''))
for (i in 2:NCOL(tmpppt_df))
{
  a = as.data.frame(bins[[i]])
  names(a) = (paste0(names(tmpppt_df)[i],'_bin',collapse = ''))
  tmpppt_df_binned = cbind (tmpppt_df_binned, a)
}



# making fia_db_with_Bins
delta= 0.1
res= 10 #on how many pieces we bin            

data_name = paste0('PRESENCE_', species_number)

fia_db[[data_name]] = (fia_db[[paste0('REL_BA_', species_number)]] > delta)

fia_db = as.data.frame(fia_db)
fia_db_with_Bins = cbind(fia_db, tmpppt_df_binned)

##################
## Data cleaning #
##################

# create a backup of the data, prior to cleaning
fia_db_with_Bins_backup = fia_db_with_Bins

# keep only the central subplot (CONDID is 1)
fia_db_with_Bins = fia_db_with_Bins[which(fia_db_with_Bins$CONDID == 1), ]

# remove duplicates, aggregate removes duplicates
form = as.formula(paste(data_name, '~ .')) # create the formula
fia_db_with_Bins = aggregate(form, fia_db_with_Bins[, c(8,9,12:ncol(fia_db_with_Bins))], max)

#nrow(fia_db_with_Bins) = 271 846

#computes intermediate variables to speed up the following computations
n_T = length(which(fia_db_with_Bins[, data_name]==T))
#n_T = 4 959

n_all = nrow(fia_db_with_Bins)
#n_all = 271 846

fia_db_with_Bins_T = fia_db_with_Bins[which(fia_db_with_Bins[, data_name]==T), ]


  
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




################################################################################################
# TEST OF THE MODEL WITH/WITHOUT INTERACTIONS, COMPUTED WITH THE EXHAUSTIVE/RESAMPLED APPROACH #
################################################################################################


# Get a subset of the climatic data for Species 
#fia_db_with_Bins = fia_db_with_Bins[sample.int(nrow(fia_db_with_Bins),50000),] # keep only a subset, for faster testing
clim_data = clim_data_generator(
  fia_db_with_Bins = fia_db_with_Bins,
  fia_db_with_Bins_T = fia_db_with_Bins[which(fia_db_with_Bins[, data_name]==T), ],
  n_all = nrow(fia_db_with_Bins),
  n_T = length(which(fia_db_with_Bins[, data_name]==T)),
  res = res,
  variables_nb = 19 # keep only a few variables, for faster testing (should be 19 for all climatic variables)
)


prediction_with = potential_with(c(1:19), clim_data)

prediction_without = potential_without(c(1:19), clim_data)

# Resampling approach to approximate the shapley values, model WITH interactions
shapleys_WITH_resampled = matrix(NA, ncol=clim_data$variables_nb, nrow=10000)
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
shapleys_WITHOUT_resampled = matrix(NA, ncol=clim_data$variables_nb, nrow=10000)
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

#prediction_with = potential_with(c(1:19), clim_data)
#prediction_without = potential_without(c(1:19), clim_data)

require('maps')
par(mfrow=c(1,3))
map('usa')
points(x=fia_db$LON[which(fia_db_with_Bins[, data_name]==T), ], 
       y=fia_db$LAT[which(fia_db_with_Bins[, data_name]==T), ])
require('graphics')
#text(-100, 40,"this is the map of USA",xpd=T)
mtext(" Eastern White Pine (Species 125) distribution", side=1)


map('usa')
points(x=fia_db$LON[which(prediction_with == T)], 
       y=fia_db$LAT[which(prediction_with == T)])
require('graphics')
#text(-100, 40,"this is the map of USA",xpd=T)
mtext("Model 'With interactions'", side=1)


map('usa')
points(x=fia_db$LON[which(prediction_without == T)], 
       y=fia_db$LAT[which(prediction_without == T)])
require('graphics')
#text(-100, 40,"this is the map of USA",xpd=T)
mtext("Model 'Without interactions'", side=1)

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