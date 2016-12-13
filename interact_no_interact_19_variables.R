code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code/' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)


require('maps')
require('raster')
require('sp')
require('rgdal')


# 125 is for Red pines (Pinus resinosa)
species_number = 121


# extraction of all the plots in the FIA DB and computation of the relative BA of the species
 source('fia_extract.R')
 fia_db = fia_extract(p=data_path, species_number)
 saveRDS(fia_db, file = paste0(code_path, '/fia_db_', species_number, '.rds'))
#fia_db = readRDS(file = paste0(code_path, '/fia_db_', species_number, '.rds'))


# Extraction of the climatic variables from Worldclim
#source('clim_extract.R')
#tmpppt = extract_tmpppt (fia_db, p=data_path) 
#saveRDS(tmpppt, file = paste0(code_path, '/tmpppt.rds'))
tmpppt = readRDS(file = paste0(code_path, '/tmpppt.rds'))



tmpppt_df = as.data.frame(tmpppt)
#summary(tmpppt_df)

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


#finding potential area for model with no interaction

fia_db_with_Bins_marker_T = rep(0, n_T)

fia_db_with_Bins_marker_all = rep(0, n_all)

dimensions = paste0('Bio', c(1:19), 's_bin')


for (i in 1:10) 
{
  fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[, dimensions[i]]*res^(i-1)
  fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[, dimensions[i]]*res^(i-1)
}

for (i in 11:19) 
{
 fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[, dimensions[i]]/res^(i-10)
 fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[, dimensions[i]]/res^(i-10)
}


#potential_area is in how many fia plots the tree can potentially grow (Interaction model)
 
potential_area = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
#potential_area = 68 623



#plot the Realized and Potential area
a = fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T))
b = fia_db_with_Bins[which(a == T), ]

par(mfrow=c(1,2))
#Realized area
map('usa')
points(x=fia_db_with_Bins_T$LON, y=fia_db_with_Bins_T$LAT)

#Potential area
map('usa')
points(x=b$LON, y=b$LAT)
 


####################################### Computing some Shapley values#####################################################################
#fia_db_with_bins has 271 846 rows
shapleys = matrix(NA, ncol=19, nrow=271846)

max_comb_size = 18  #the max size of the group of variables (the group is different from the set of all 19 variables)

pb = txtProgressBar(max=nrow(shapleys), style = 3)

for (trial in 1:nrow(shapleys)) {
  
  setTxtProgressBar(pb, trial)
  
  for (var in 1:ncol(shapleys)) {
    
    n_vars = sample.int(max_comb_size, 1) #the size of selected group  ?chosen randomly??? there can be repetitions
    
    vars_selected = c(sample((1:19)[-var], n_vars), var)
    
    dimensions = paste0('Bio', vars_selected, 's_bin')
    
    names(dimensions) = paste0('dim', 1:(n_vars+1))
    
    fia_db_with_Bins_marker_T = rep(0, n_T)
    
    fia_db_with_Bins_marker_all = rep(0, n_all)
    
    for (i in 1:n_vars) {
      
      if (i <= 10) {
        
        # working around precision/max value limitation
        
        fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]*res^(i-1)
        
        fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]*res^(i-1)
        
      } else {
        
        fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]/res^(i-10)
        
        fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]/res^(i-10)
        
      }
      
    }
    
    score_without = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
    
    # compute the score with the additional var:
    
    i = i + 1
    
    if (i <= 10) {
      
      # working around precision/max value limitation
      
      fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]*res^(i-1)
      
      fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]*res^(i-1)
      
    } else {
      
      fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]/res^(i-10)
      
      fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]/res^(i-10)
      
    }
    
    score_with = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
#score_without = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
    
    shapleys[trial, var] = (score_without - score_with) * ( factorial(n_vars) * factorial(max_comb_size-n_vars) / factorial(19) )
    
  }
  
}


saveRDS(shapleys, file = paste0(code_path, '/shapleys_', species_number, '.rds'))

n_done = which(is.na(shapleys[,ncol(shapleys)]))[1] - 1
#n_done = 22 754

###################################another approach to compute Shapley Values#########################################################
#2^19-1 = 524 287
shap_matrix = matrix(NA, ncol=19, nrow = 524287)

for (k in 1:524287) 
{ 
  #print(intToBits(k)) 
  vars_selected=c()
  for (j in 1:19) 
  {
    if (intToBits(k)[j] == 01) 
    {
      vars_selected = c(vars_selected, j) 
    }
    
  }
  #print(vars_selected)
  vars_selected_complement = (1:19)[-vars_selected]
  #print(vars_selected_complement)
  
  for (var in 1:19) 
   {
     if (var %in% vars_selected == FALSE) 
      { source('Shapleys_function.R')
        shapleys = shapleys(vars_selected, var)
        shap_matrix[k, var]= shapleys }
    
    else {vars_selected = vars_selected_complement
          source('Shapleys_function.R')
          shapleys = shapleys(vars_selected, var)
          shap_matrix[k, var]= shapleys }
    } 
}


saveRDS(shap_matrix, file = paste0(code_path, '/shap_matrix_', species_number, '.rds'))

n_done = which(is.na(shap_matrix[,ncol(shap_matrix)]))[1] - 1
#n_done = 16 139



#################### Code to output Shapley Values  #####################################
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
                       
                       'Precipitation of Coldest Quarter')


bioclim_types = c(1, 3, 3, 3, 1, 1, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2) # this variable encodes whether the bioclimatic variable is linked to temperature (1), precipitation (2) or if it is some measure of variability (3)

colormap = c('orange','lightblue', 'purple')


species_name= 'Sitka Spruce'


## THE EXACT SHAPLEY VALUES:

shapleys = readRDS('C:/Users/Olga Rumyantseva/Documents/R Files/Code/shap_matrix_125.rds')


scores_exact = colMeans(shapleys, na.rm=T)

order_exact = order(scores_exact, decreasing = T)

bp = barplot(scores_exact[order_exact], col=colormap[bioclim_types[order_exact]], names.arg=NA, las=1, ylab='Exact shapley value')

text(bp, 1, all_variable_names[order_exact], srt=90, adj=0, xpd=T, font=3)

legend('topright', fill=colormap, c('Temperature-related', 'Precipitation-related', 'Variability-related'), bty='n', xpd=T)

title(paste0('Which bioclimatic variables characterize best the potential area of ', species_name, '?'))



## Approximated shapley values (with random sampling):

shapleys = readRDS('C:/Users/Olga Rumyantseva/Documents/R Files/Code/shap_matrix_125.rds')

n_done = which(is.na(shapleys[,ncol(shapleys)]))[1] - 1

#half_done = n_done/2


#mean1 = colMeans(shapleys[1:half_done,])

#mean2 = colMeans(shapleys[(1+half_done):n_done,])

#col_diff = mean1/mean2*100-100

#total_diff = mean(col_diff)

#abs_col_diff = abs(col_diff)



# names(col_diff) = c('Annual Mean Temperature',
#                    
#                    'Mean Diurnal Range (Mean of monthly (max temp - min temp))',
#                    
#                    'Isothermality (BIO2/BIO7) (* 100)',
#                    
#                    'Temperature Seasonality (standard deviation *100)',
#                    
#                    'Max Temperature of Warmest Month',
#                    
#                    'Min Temperature of Coldest Month',
#                    
#                    'Temperature Annual Range (BIO5-BIO6)',
#                    
#                    'Mean Temperature of Wettest Quarter',
#                    
#                    'Mean Temperature of Driest Quarter',
#                    
#                    'Mean Temperature of Warmest Quarter',
#                    
#                    'Mean Temperature of Coldest Quarter',
#                    
#                    'Annual Precipitation',
#                    
#                    'Precipitation of Wettest Month',
#                    
#                    'Precipitation of Driest Month',
#                    
#                    'Precipitation Seasonality (Coefficient of Variation)',
#                    
#                    'Precipitation of Wettest Quarter',
#                    
#                    'Precipitation of Driest Quarter',
#                    
#                    'Precipitation of Warmest Quarter',
#                    
#                    'Precipitation of Coldest Quarter')




###### names(abs_col_diff) = c('Annual Mean Temperature',
#                        
#                        'Mean Diurnal Range (Mean of monthly (max temp - min temp))',
#                        
#                        'Isothermality (BIO2/BIO7) (* 100)',
#                        
#                        'Temperature Seasonality (standard deviation *100)',
#                        
#                        'Max Temperature of Warmest Month',
#                        
#                        'Min Temperature of Coldest Month',
#                        
#                        'Temperature Annual Range (BIO5-BIO6)',
#                        
#                        'Mean Temperature of Wettest Quarter',
#                        
#                        'Mean Temperature of Driest Quarter',
#                        
#                        'Mean Temperature of Warmest Quarter',
#                        
#                        'Mean Temperature of Coldest Quarter',
#                        
#                        'Annual Precipitation',
#                        
#                        'Precipitation of Wettest Month',
#                        
#                        'Precipitation of Driest Month',
#                        
#                        'Precipitation Seasonality (Coefficient of Variation)',
#                        
#                        'Precipitation of Wettest Quarter',
#                        
#                        'Precipitation of Driest Quarter',
#                        
#                        'Precipitation of Warmest Quarter',
#                        
#                        'Precipitation of Coldest Quarter')




require(boot) ## If error, install: install.packages("boot")

get_mean_subsample = function(data, indices) {
  
  return(colMeans(data[indices, ]))
  
}

results = boot(data=shapleys[1:n_done, ], statistic=get_mean_subsample, R=100)

# plot(results)

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

bioclim_types = c(1, 3, 3, 3, 1, 1, 3, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2) # this variable encodes whether the bioclimatic variable is linked to temperature (1), precipitation (2) or if it is some measure of variability (3)

ordered_shap = confs[, order(confs[1,], decreasing = T)]

scores = unlist(ordered_shap[1,])

colormap = c('orange','lightblue', 'purple')




# plot the figure

bp = barplot(scores, col=colormap[bioclim_types[unlist(ordered_shap[5,])]], names.arg=NA, las=1, ylab='Approximated shapley value')

text(bp, 1, names(ordered_shap[1,]), srt=90, adj=0, xpd=T, font=3)

legend('topright', fill=colormap, c('Temperature-related', 'Precipitation-related', 'Variability-related'), bty='n', xpd=T)

title(paste0('Which bioclimatic variables characterize best the potential area of ', species_name, '?'))


require('Hmisc')

errbar(bp, unlist(ordered_shap[1,]), unlist(ordered_shap[3,]), unlist(ordered_shap[4,]), add=T, xpd=T, pch=NA)









