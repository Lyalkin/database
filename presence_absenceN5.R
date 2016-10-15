#"Main with N" placeholder


code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('maps')
require('raster')
require('sp')


BIO1 = T ; names(BIO1) = "Annual Mean Temperature"
BIO2 = T ; names(BIO2) = "Mean Diurnal Range (Mean of monthly (max temp - min temp))"
BIO3 = T ; names(BIO3) = "Isothermality (BIO2/BIO7) (* 100)"
BIO4 = T ; names(BIO4) = "Temperature Seasonality (standard deviation *100)"
BIO5 = T ; names(BIO5) = "Max Temperature of Warmest Month"
BIO6 = T ; names(BIO6) = "Min Temperature of Coldest Month"
BIO7 = T ; names(BIO7) = "Temperature Annual Range (BIO5-BIO6)"
BIO8 = T ; names(BIO8) = "Mean Temperature of Wettest Quarter"
BIO9 = T ; names(BIO9) = "Mean Temperature of Driest Quarter"
BIO10 = T ; names(BIO10) = "Mean Temperature of Warmest Quarter"
BIO11 = T ; names(BIO11) = "Mean Temperature of Coldest Quarter"
BIO12 = T ; names(BIO12) = "Annual Precipitation"
BIO13 = T ; names(BIO13) = "Precipitation of Wettest Month"
BIO14 = T ; names(BIO14) = "Precipitation of Driest Month"
BIO15 = T ; names(BIO15) = "Precipitation Seasonality (Coefficient of Variation)"
BIO16 = T ; names(BIO16) = "Precipitation of Wettest Quarter"
BIO17 = T ; names(BIO17) = "Precipitation of Driest Quarter"
BIO18 = T ; names(BIO18)= "Precipitation of Warmest Quarter"
BIO19 = T ; names(BIO19) = "Precipitation of Coldest Quarter"

clim_vars_in_use = c(BIO1,BIO2,BIO3,BIO4,BIO5,BIO6,BIO7,BIO8,BIO9,
                     BIO10,BIO11,BIO12,BIO13,BIO14,BIO15,BIO16,BIO17,BIO18,BIO19)
thresh = 0.1
res = 10
trials=1
tmpppt = readRDS(file='tmpppt_98.rds') #always 98 because it's only needed once built in bingeneration.R
species_list = read.csv('ref_FIA_QB_updated_withcurtis.csv')
species_list = species_list[order(species_list[,3]),]

###########################################################
####################set CCV here###########################

data_number = 212    #was data_number  98

###########################################################
###########################################################
# 
 #source('fia_extract.R')
 #fia_db = fia_extract(data_path, species_to_report=data_number)
 #saveRDS(fia_db,file=paste0('fia_db_',data_number,'.rds'))

######################
######################

# Jean commented these lines on 9-1-2016
# row_number = which(species_list[,3] == data_number, arr.ind=TRUE)
# species_name = droplevels(species_list[[row_number,5]])
data_name = paste0('PRESENCE_',data_number)
fia_db = readRDS(file = paste0(code_path, '/fia_db_', species_number, '.rds'))  #olga added
#fia_db = readRDS(file=paste0('fia_db_',data_number,'.rds'))
fia_db[[data_name]] = (fia_db[[paste0('REL_BA_',data_number)]] > thresh)
Plot_count = sum(fia_db[[data_name]],na.rm = T)
map('state')
points(fia_db$LON, fia_db$LAT, cex=fia_db[[data_name]])

tmppptDF = as.data.frame(tmpppt) # replaces the fia_dp$tmp = tmpppt$tmps... etc.
fia_dbN = cbind(fia_db, tmppptDF)

source("BinGeneration.R")
tmpppt_Bins = BinGeneration(tmppptDF,res)
fia_db_with_Bins = cbind(fia_db,tmpppt_Bins)


# cropping the worldwide rasters to the extent of the US

load('mask_usa.RData')
extent_limit = extent(mask_usa)
extents_name= "USA"



##########################
# 
# source("BinGeneration.R")
# tmpppt_raster_bins = BinGeneration_Raster(tmppptDF,res,extent_limit,extents_name,clim_vars_in_use)
# saveRDS(tmpppt_raster_bins, file=paste0('tmpppt_raster_bins_',data_number,'.rds'))

##########################


tmpppt_raster_bins = readRDS(file=paste0('tmpppt_raster_bins_',data_number,'.rds'))

###################################################################################################################

######################
### Shapely Method ###
######################


# computes intermediate variables to speed up the following computations
fia_db_with_Bins_trimmed = aggregate(as.formula(paste(data_name, '~ .')), fia_db_with_Bins[, c(8,9,12:ncol(fia_db_with_Bins))], max)
# The line above looks for repeated surveys of the same plot (with the same LAT, LON)
# it keeps only one line per survey
# if at least one of the repeated surveys contained the species of interest, then it keeps this survey (fia_db[[data_name]] is T)
# if none of them contained the species of interest, then it keeps only one line for the repeated measurements, with fia_db[[data_name]] being F
# advice: write this in another way, using a for loop?

# saveRDS(fia_db_with_Bins_trimmed, file=paste0('fia_db_with_Bins_trimmed_',data_number,'.rds'))
n_T = length(which(fia_db_with_Bins_trimmed[, data_name]==T))
n_all = nrow(fia_db_with_Bins_trimmed)
fia_db_with_Bins_T = fia_db_with_Bins_trimmed[which(fia_db_with_Bins_trimmed[, data_name]==T), ]



##################################


#saveRDS(shapleys, file=paste0('Shapleys_',data_number,'.rds'))
shapleys=readRDS(file=paste0('Shapleys_',data_number,'.rds'))



###########################
## Mass Species Analysis ##
###########################


completed_species = read.csv('already_completed_species.csv')

for (i in 1:14) {
  data_number= as.numeric(completed_species[1,i])
  shapleys=readRDS(file=paste0('shapleys_',data_number,'.rds'))
  row_number = which(species_list[,3] == data_number, arr.ind=TRUE)
  species_name = droplevels(species_list[[row_number,5]])
  data_name = paste0('PRESENCE_',data_number)
  
  n_done = which(is.na(shapleys[,ncol(shapleys)]))[1] - 1
  half_done = n_done/2
  mean1 = colMeans(shapleys[1:half_done,])
  mean2 = colMeans(shapleys[(1+half_done):n_done,])
  col_diff = mean1/mean2*100-100
  total_diff = mean(col_diff)
  abs_col_diff = abs(col_diff)
  names(col_diff) = c('Annual Mean Temperature',
                      'Mean Diurnal Range (Mean of monthly (max temp - min temp))',
                      'Isothermality (BIO2/BIO7) (* 100)',
                      'Temperature Seasonality (standard deviation *100)',
                      'Max Temperature of Warmest Month',
                      'Min Temperature of Coldest Month',
                      'Temperature Annual Range (BIO5-BIO6)',
                      'Mean Temperature of Wettest Quarter',
                      'Mean Temperature of Driest Quarter',
                      'Mean Temperature of Warmest Quarter',
                      'Mean Temperature of Coldest Quarter',
                      'Annual Precipitation',
                      'Precipitation of Wettest Month',
                      'Precipitation of Driest Month',
                      'Precipitation Seasonality (Coefficient of Variation)',
                      'Precipitation of Wettest Quarter',
                      'Precipitation of Driest Quarter',
                      'Precipitation of Warmest Quarter',
                      'Precipitation of Coldest Quarter')
  names(abs_col_diff) = c('Annual Mean Temperature',
                          'Mean Diurnal Range (Mean of monthly (max temp - min temp))',
                          'Isothermality (BIO2/BIO7) (* 100)',
                          'Temperature Seasonality (standard deviation *100)',
                          'Max Temperature of Warmest Month',
                          'Min Temperature of Coldest Month',
                          'Temperature Annual Range (BIO5-BIO6)',
                          'Mean Temperature of Wettest Quarter',
                          'Mean Temperature of Driest Quarter',
                          'Mean Temperature of Warmest Quarter',
                          'Mean Temperature of Coldest Quarter',
                          'Annual Precipitation',
                          'Precipitation of Wettest Month',
                          'Precipitation of Driest Month',
                          'Precipitation Seasonality (Coefficient of Variation)',
                          'Precipitation of Wettest Quarter',
                          'Precipitation of Driest Quarter',
                          'Precipitation of Warmest Quarter',
                          'Precipitation of Coldest Quarter')
  
  
  
  
  # bootstrapping with 1000 replications 
  results = boot(data=shapleys[1:n_done, ], statistic=get_mean_subsample, R=10)
  
  # Bootstrap 95% CI for R-Squared
  require(boot)
  # function to obtain the mean from the data 
  get_mean_subsample = function(data, indices) {
    return(colMeans(data[indices, ]))
  }
  # return(mean(data[indices, 1]))
  # view results
  results
  filename=paste0(data_number,'_histogram','.png')
  plot(results)
  # dev.off()
  
  # get 95% confidence interval 
  confs = sapply(1:19, function(i){
    tmp=boot.ci(results, type="norm", index=i)
    return(c(tmp$t0, tmp$normal))}
  )
  confs = rbind(confs, 1:19)
  row.names(confs) = c('statistic', 'conf_int', 'min_bound', 'max_bound', 'climatic_var_i')
  colnames(confs) = c('Annual Mean Temperature',
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
                      'Precipitation of Coldest Quarter')
  
  # visualization of the results with a colored barplot with 95% CI
  bioclim_types = c(1, 3, 3, 3, 1, 1, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2) # this variable encodes whether the bioclimatic variable is linked to temperature (1), precipitation (2) or if it is some measure of variability (3)
  ordered_shap = confs[, order(confs[1,], decreasing = T)]
  
  scores = unlist(ordered_shap[1,])
  
  colormap = c('orange','lightblue', 'purple')
  png(filename=paste0(data_number,'_data_rankings','.png'))
  bp = barplot(scores, col=colormap[bioclim_types[unlist(ordered_shap[5,])]], names.arg=NA, las=1, ylab='Shapley value')
  text(bp, 1, names(ordered_shap[1,]), srt=90, adj=0, xpd=T, font=3)
  legend('topright', fill=colormap, c('Temperature-related', 'Precipitation-related', 'Variability-related'), bty='n', xpd=T)
  title(paste0('Which bioclimatic variables characterize best the potential area of ', species_name, '?'))
  require('Hmisc')
  errbar(bp, unlist(ordered_shap[1,]), unlist(ordered_shap[3,]), unlist(ordered_shap[4,]), add=T, xpd=T, pch=NA)
  dev.off()
  ## correlation matrix (might take a long time for widespread species)
  require('corrgram')
  png(filename=paste0(data_number,'_corrgram','.png'))
  corrgram(tmppptDF[fia_db[[data_name]], unlist(ordered_shap[5,])], upper.panel = panel.ellipse, lower.panel = panel.conf, labels = colnames(ordered_shap), cex.labels = 0.5)
  dev.off()
  
  
  #######
  # PCA #
  #######
  
  get_PCA_vars = function(tmppptDF) {
    not_na = which(complete.cases(tmppptDF))
    pca = prcomp(tmppptDF[not_na, ], center = T, scale. = T)
    pcaDF = matrix(NA, nrow=nrow(tmppptDF), ncol=ncol(tmppptDF))
    pcaDF[not_na, ] = pca$x
    pcaDF = data.frame(pcaDF)
    colnames(pcaDF) = colnames(tmppptDF)
    return(pcaDF)
  }
  
  pcaDF = get_PCA_vars(tmppptDF)
  
  
  ## PCA QAD: on the whole dataset
  pca = prcomp(tmppptDF[complete.cases(tmppptDF), ], center = T, scale. = T)
  ## PCA QAD: only for Eastern Hemlocks
  where_species_live = which((fia_dbN[[data_name]]==1) & (complete.cases(tmppptDF)))
  pca = prcomp(tmppptDF[where_species_live, ], center = T, scale. = T)
  
  pca_trimmed = pca
  pca_trimmed$x = pca_trimmed$x[sample.int(nrow(pca_trimmed$x), 1000),]
  
  pdf(file= paste0(species_name,'_PCA_biplot'), height = 9, width = 9)
  biplot(pca_trimmed)
  title('biplot of a subsample of 1000 data points in PC1/PC2 space')
  dev.off()
  
  pdf(file=paste0(species_name,'_PCA_Cumulative_Variance'), height = 9, width = 9)
  barplot(cumsum(pca_trimmed$sdev) / sum(pca_trimmed$sdev))
  title('Cumulative proportion of variance explained by each Principal Components')
  abline(h=0.9)
  dev.off()
}







# 
# get_PCA_vars = function(tmppptDF) {
#   not_na = which(complete.cases(tmppptDF))
#   pca = prcomp(tmppptDF[not_na, ], center = T, scale. = T)
#   pcaDF = matrix(NA, nrow=nrow(tmppptDF), ncol=ncol(tmppptDF))
#   pcaDF[not_na, ] = pca$x
#   pcaDF = data.frame(pcaDF)
#   colnames(pcaDF) = colnames(tmppptDF)
#   return(pcaDF)
# }
# 
# pcaDF = get_PCA_vars(tmppptDF)
# 
# 
# ## PCA QAD: on the whole dataset
# pca = prcomp(tmppptDF[complete.cases(tmppptDF), ], center = T, scale. = T)
# ## PCA QAD: only for Eastern Hemlocks
# where_species_live = which((fia_dbN[[data_name]]==1) & (complete.cases(tmppptDF)))
# pca = prcomp(tmppptDF[where_species_live, ], center = T, scale. = T)
# 
# pca_trimmed = pca
# pca_trimmed$x = pca_trimmed$x[sample.int(nrow(pca_trimmed$x), 1000),]
# 
# biplot(pca_trimmed)
# title('biplot of a subsample of 1000 data points in PC1/PC2 space')
# 
# # barplot(pca_trimmed$sdev / sum(pca_trimmed$sdev))
# # title('Proportion of variance explained by each Principal Components')
# 
# barplot(cumsum(pca_trimmed$sdev) / sum(pca_trimmed$sdev))
# title('Cumulative proportion of variance explained by each Principal Components')
# abline(h=0.9)
# 
