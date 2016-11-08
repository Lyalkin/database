code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

#install.packages("scatterplot3d", dependencies = TRUE)
#install.packages("rgl", dependencies = TRUE)

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

#computes intermediate variables to speed up the following computations
n_T = length(which(fia_db_with_Bins[, data_name]==T))
#n_T = 12 414

n_all = nrow(fia_db_with_Bins)
#n_all = 754 851

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


#potential_area is in how many fia plots the tree can potentially grow
 
potential_area = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
#potential_area = 173 566



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
#fia_db   has 754 851 rows
#shapleys has 1 000 000 rows
shapleys = matrix(NA, ncol=19, nrow=1000000)

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







