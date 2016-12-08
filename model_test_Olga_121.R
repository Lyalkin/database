# set the working directory   lines 2, 4, 5, 16, 135, 206
setwd("~/tmptodelete/model_without_interactions/") # change this on your computer

code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code/' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)


require('maps')
require('raster')
require('sp')
require('rgdal')


# Pinus Palustris,  Longleaf Pine
species_number = 121


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


#delta = 0
res = 10 #on how many pieces we bin            

#pl = which(fia_db$REL_BA_121 > 0)
#require('maps')
#map('usa')
#points(x=fia_db$LON[pl], y=fia_db$LAT[pl])

data_name = paste0('PRESENCE_', species_number)

fia_db[[data_name]] = (fia_db[[paste0('REL_BA_', species_number)]] > 0)

fia_db = as.data.frame(fia_db)

fia_db_with_Bins = cbind(fia_db, tmpppt_df_binned)
#head(fia_db_with_Bins)

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
#q are plots where our tree present
q = which(fia_db_with_Bins[, data_name]==T)
n_T = length(q)
n_T 

n_all = nrow(fia_db_with_Bins)
n_all 

fia_db_with_Bins_T = fia_db_with_Bins[q, ]

  
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

#q = which(fia_db_with_Bins[, data_name]==T)
q_int  = which(prediction_with==T)
q_nint = which(prediction_without==T)

require('maps')
par(mfrow=c(1,3))
map('usa')
points(x=fia_db_with_Bins$LON[q], y=fia_db_with_Bins$LAT[q])
require('graphics')
mtext(" Longleaf Pine (Pinus Palustris 121) distribution", side=1)

map('usa')
points(x=fia_db_with_Bins$LON[q_int], y=fia_db_with_Bins$LAT[q_int])
require('graphics')
mtext("Model 'With interactions'", side=1)


map('usa')
points(x=fia_db_with_Bins$LON[q_nint], y=fia_db_with_Bins$LAT[q_nint])
require('graphics')
mtext("Model 'Without interactions'", side=1)

