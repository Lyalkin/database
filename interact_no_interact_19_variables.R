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
summary(tmpppt_df)

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
data_name = paste0('PRESENCE_', species_number)

fia_db[[data_name]] = (fia_db[[paste0('REL_BA_', species_number)]] > 0.1)

fia_db = as.data.frame(fia_db)
fia_db_with_Bins = cbind(fia_db, tmpppt_df_binned)


# computes intermediate variables to speed up the following computations
n_T = length(which(fia_db_with_Bins[, data_name]==T))
#n_T = 12 414

n_all = nrow(fia_db_with_Bins)
#n_all = 754 851

fia_db_with_Bins_T = fia_db_with_Bins[which(fia_db_with_Bins[, data_name]==T), ]


# Computing some Shapley values
#fia_db has 754 851 rows
#        1 000 000
shapleys = matrix(NA, ncol=19, nrow=1000000)

max_comb_size = 18

pb = txtProgressBar(max=nrow(shapleys), style = 3)

for (trial in 1:nrow(shapleys)) {
  
  setTxtProgressBar(pb, trial)
  
  for (var in 1:ncol(shapleys)) {
    
    n_vars = sample.int(max_comb_size, 1)
    
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
    
    shapleys[trial, var] = (score_without - score_with) * ( factorial(n_vars) * factorial(max_comb_size-n_vars) / factorial(19) )
    
  }
  
}



plots_with_tree = which(fia_db[[paste0("REL_BA_", species_number)]] > 0)


# p is a presence_absence vector for FIA plots if we consider model with no interactions

d = dim(tmpppt_df_binned)[1]
p = vector(mode = "logical", length = d)  #d= -number of fia plots, vector has FALSE components initially

c = as.data.frame(cbind(tmpppt_df_binned[,var1], 
                        tmpppt_df_binned[,var2], tmpppt_df_binned[,var3]))


for (i in 1:d) {p[i] = vars_no_interact[c[i,1], c[i,2], c[i,3]]}

for (i in 1:d) 
{
  if (is.na(p[i]) == TRUE) {p[i] = 0} 
}

# p2 is a presence_absence vector for FIA plots if we consider model with interactions

p2 = vector(mode = "logical", length = d) # prediction of the model with interaction

for (i in 1:d) {
  
  p2[i] = vars[c[i,1], c[i,2], c[i,3]]
  
  if (is.na(p2[i]) == TRUE) {
    
    p2[i] = 0
    
  }
}

potential_area_with_interactions = sum(p2)
#potential_area_with_interactions = 383 950 

# finding Potential area for the tree (in how many fia plots the tree can be, if we consider
#the model without interactions)

potential_area_no_interactions = sum(p)
#potential_area_no_interactions = 384 979

# finding Relative Basal Area for 
#plots_with_tree = which(fia_db[[paste0("REL_BA_", species_number)]] > 0)

total_basal_area = sum(fia_db$BASAL.AREA[plots_with_tree]) 
#total_basal_area =  22 406.7




#plot together Rel Bas Area, Potential Area(no ineract model), Potential Area(ineract model) for Giant Seq
p = as.data.frame(p)
p2 = as.data.frame(p2)

plots_no_interact = which(p > 0)
plots_interact = which(p2 > 0)




