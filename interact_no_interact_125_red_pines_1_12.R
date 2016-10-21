code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('maps')
require('raster')
require('sp')
require('rgdal')

# 125 is for Red pines (Pinus resinosa)
species_number = 125

#BIO5 = "Max Temperature of Warmest Month"
#BIO6 = "Min Temperature of Coldest Month"
var1 = 1      #climatic variable 1
var2 = 12     #climatic variable 2

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


#BIO1= Annual Mean Temperature
#BIO12= Annual Precipitation

plots_with_species = which(fia_db[paste0('REL_BA_', species_number)] > 0)


# Now cut (make bins) for all NCOL(tmpppt_df) = 19 climatic variables
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
#head(tmpppt_df_binned)
#dim(tmpppt_df_binned)




# we can plot clim variables


plot(tmpppt_df_binned[,var1][plots_with_tree], tmpppt_df_binned[,var2][plots_with_tree], type = "p",
     xlab=paste0("clim var ", var1), ylab=paste0("clim var ", var2), 
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)




#var1_var2_tree has columns var1 var2 binned for which the given tree grows

plots_with_tree = which(fia_db[[paste0("REL_BA_", species_number)]] > 0)

var1_var2_tree = cbind(tmpppt_df_binned[,var1][plots_with_tree], tmpppt_df_binned[,var2][plots_with_tree])
var1_var2_tree = na.omit(unique(as.data.frame(var1_var2_tree)))
#var1_var2_tree


#vars is a matrix having vars[i,j]=1 if var1_var2_tree has a row "i j" 
#and                     vars[i,j]=0 otherwise

vars = matrix(0, 10, 10)
for (i in 1:dim(var1_var2_tree)[1])
{vars[var1_var2_tree[i,1],var1_var2_tree[i,2]] = 1}
#vars



#vars_no_interact obtained from vars matrix making it 'looking convex'
#a1, a2 -  matrices for loop working

vars_no_interact = vars
vars_no_interact 

a1= matrix(7, 10, 10)
a2 = matrix(77, 10, 10)


while (identical(a1, a2) == FALSE) 
{  
  a1 = vars_no_interact
  print(a1)
  
  for (i in 1:9)
  {
    for (j in 1:9)  
    {
      if (vars_no_interact[i,j] ==   0 & vars_no_interact[i,j+1] ==   1 & 
          vars_no_interact[i+1,j] == 1 & vars_no_interact[i+1,j+1] == 1) 
        
      {vars_no_interact[i,j] = 1}
      
      
      if  (vars_no_interact[i,j] ==   1 & vars_no_interact[i,j+1] ==   0 & 
           vars_no_interact[i+1,j] == 1 & vars_no_interact[i+1,j+1] == 1) 
        
      {vars_no_interact[i,j+1] = 1}
      
      
      if (vars_no_interact[i,j] ==   1 & vars_no_interact[i,j+1]   == 1 & 
          vars_no_interact[i+1,j] == 0 & vars_no_interact[i+1,j+1] == 1) 
        
      {vars_no_interact[i+1,j] = 1}
      
      if (vars_no_interact[i,j] ==   1   & vars_no_interact[i,j+1]   == 1 & 
          vars_no_interact[i+1,j] == 1 & vars_no_interact[i+1,j+1] ==   0) 
        
      {vars_no_interact[i+1,j+1] = 1}
    }
  }
  a2 = vars_no_interact
  print(a2)
}  

print(vars_no_interact)



# p is a presence_absence vector for FIA plots if we consider model with no interactions

d = dim(tmpppt_df_binned)[1]
p = vector(mode = "logical", length = d)  #d= -number of fia plots, vector has FALSE components initially

c = as.data.frame(cbind(tmpppt_df_binned[,var1], tmpppt_df_binned[,var2]))

for (i in 1:d) {p[i] = vars_no_interact[c[i,1], c[i,2]]}

for (i in 1:d) 
{
  if (is.na(p[i]) == TRUE) {p[i] = 0} 
}

# p2 is a presence_absence vector for FIA plots if we consider model with interactions

p2 = vector(mode = "logical", length = d) # prediction of the model with interaction

for (i in 1:d) {
  
  p2[i] = vars[c[i,1], c[i,2]]
  
  if (is.na(p2[i]) == TRUE) {
    
    p2[i] = 0
    
  }
}

potential_area_with_interactions = sum(p2)
#potential_area_with_interactions = 444 185

# finding Potential area for the tree (in how many fia plots the tree can be, if we consider
#the model without interactions)

potential_area_no_interactions = sum(p)
#potential_area_no_interactions = 504 272

# finding Relative Basal Area for 
#plots_with_tree = which(fia_db[[paste0("REL_BA_", species_number)]] > 0)

total_basal_area = sum(fia_db$BASAL.AREA[plots_with_tree]) 
#total_basal_area = 22 406.7




#plot together Rel Bas Area, Potential Area(no ineract model), Potential Area(ineract model) for Giant Seq
p = as.data.frame(p)
p2 = as.data.frame(p2)

plots_no_interact = which(p > 0)
plots_interact = which(p2 > 0)

require('maps')

par(mfrow=c(1,3))
map('usa')
points(x=fia_db$LON[plots_with_tree], y=fia_db$LAT[plots_with_tree])

map('usa')
points(x=fia_db$LON[plots_no_interact], y=fia_db$LAT[plots_no_interact])

map('usa')
points(x=fia_db$LON[plots_interact], y=fia_db$LAT[plots_interact])


###########################  some testing of the code  ###########################################

cat("Annual temperature / 
    precipitation of 30 places predicted to be giant-sequoia-friendly (potential area)\n")

random_subset_30 = sample(which(p>0), size=30)

for (random_place in random_subset_30) {
  
  cat("temp: ", tmpppt_df[random_place, 1], " & ppt: ", tmpppt_df[random_place, 12], "\n")
  
}


cat("Annual temperature / precipitation of places where we know 
    giant sequoias are there (realized area)\n")


for (place_with_tree in which(fia_db$REL_BA_species_number > 0)) {
  
  cat("temp: ", tmpppt_df[place_with_tree, 1], " & ppt: ", tmpppt_df[place_with_tree, 12], "\n")
  
}

#-----------------------------------------------binned----------------------------------------------


for (random_place in random_subset_30) {
  
  cat("temp: ", tmpppt_df_binned[random_place, 1], " & ppt: ", tmpppt_df_binned[random_place, 12], "\n")
  
}


cat("Annual temperature / precipitation of places where we know 
    giant sequoias are there (realized area)\n")


for (place_with_tree in which(fia_db$REL_BA_species_number > 0)) {
  
  cat("temp: ", tmpppt_df_binned[place_with_tree, 1], " & ppt: ", tmpppt_df_binned[place_with_tree, 12], "\n")
  
}

#--------------------------------------------------------------------------------------------
s = which(tmpppt_df_binned[,1] == 5 & tmpppt_df_binned[,12] == 4)
v = as.data.frame(tmpppt_df_binned[s, ])
dim(v)
#---------------------------------------------------------------------------------------------



