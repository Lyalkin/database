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

#BIO5 = "Max Temperature of Warmest Month"
#BIO6 = "Min Temperature of Coldest Month"
#BIO7 = "Temperature Annual Range (BIO5-BIO6)"
var1 = 5      #climatic variable 1
var2 = 6      #climatic variable 2
var3 = 7      #climatic variable 3

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



plots_with_tree = which(fia_db[paste0('REL_BA_', species_number)] > 0)


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


# 3D plots
library(rgl)

x1 = tmpppt_df[,var1][plots_with_tree] 
#summary(x1)
y1 = tmpppt_df[,var2][plots_with_tree] 
#summary(y1)
z1 = tmpppt_df[,var3][plots_with_tree]   
#summary(z1)

x = tmpppt_df_binned[,var1][plots_with_tree]
y = tmpppt_df_binned[,var2][plots_with_tree]
z = tmpppt_df_binned[,var3][plots_with_tree]


plot3d(tmpppt_df_binned[,var1][plots_with_tree], tmpppt_df_binned[,var2][plots_with_tree],
       tmpppt_df_binned[,var3][plots_with_tree], 
       xlim=c(0, 11), ylim=c(0, 11), zlim=c(0, 11),
       xlab=paste0("clim var ", var1), ylab=paste0("clim var ", var2), zlab=paste0("clim var ", var3), col="red", size=3)


plot3d(x1, y1, z1, 
              xlim=c(220, 335), ylim=c(-250, -50), zlim=c(300, 505),
              xlab=paste0("clim var ", var1), ylab=paste0("clim var ", var2), zlab=paste0("clim var ", var3),
              highlight.3d=TRUE, col.axis="blue",
              col.grid="lightblue",  pch=20)




#v has columns var1 var2 var3 binned for which the given tree grows

plots_with_tree = which(fia_db[[paste0("REL_BA_", species_number)]] > 0)

v = cbind(tmpppt_df_binned[,var1][plots_with_tree], 
          tmpppt_df_binned[,var2][plots_with_tree],
          tmpppt_df_binned[,var3][plots_with_tree])
v = na.omit(unique(as.data.frame(var1_var2_var3_tree)))



#dim(v) = 30 rows 3 columns

vars = array(0, dim=c(10,10,10))
vars_no_interact = array(0, dim=c(10,10,10))

v1 = unique(v[ ,1])
v2 = unique(v[ ,2])
v3 = unique(v[ ,3])



source('matrix_inter_no_inter.R')

#m has columns v1 v2 v3 where the column v1 has the same values  
for (k in 1:length(v1))
{ 
 n = which(v[ ,1] == v1[k])
 m = unique(v[n, ])
 print(m)
  for (i in 1:dim(m)[1])  {vars[v1[k], m[i,2], m[i, 3]] = 1}
  
  print(paste0('matrix with interactions vars[', v1[k],  ', , ]'))
  print(vars[v1[k], , ])
 
  vars_no_interact[v1[k], , ] = conv(vars[v1[k], , ])
  
  print(paste0('matrix without interactions vars_no_interact[', v1[k],  ', , ]'))
  print(vars_no_interact[v1[k], , ])
}

#m has columns v1 v2 v3 where the column v2 has the same values
for (k in 1:length(v2))
{ 
  n = which(v[ ,2] == v2[k])
  m = unique(v[n, ])
  print(m)
  
  for (i in 1:dim(m)[1])  {vars[m[i,1], v2[k], m[i, 3]] = 1}
  
  print(paste0('matrix with interactions vars[  ,', v2[k],', ]'))
  print(vars[ ,v2[k], ])
  
  vars_no_interact[ ,v2[k], ] = conv(vars[ ,v2[k], ])
  
  print(paste0('matrix with interactions vars_no_interact[  ,', v2[k],', ]'))
  print(vars_no_interact[ ,v2[k], ])
}

#m has columns v1 v2 v3 where the column v3 has the same values
for (k in 1:length(v3))
{ 
  n = which(v[ ,3] == v3[k])
  m = unique(v[n, ])
  print(m)
  
  for (i in 1:dim(m)[1])  {vars[m[i,1], m[i, 3], v3[k]] = 1}
  
  print(paste0('matrix with interactions vars[  , ,', v3[k],']'))
  print(vars[ , ,v3[k]])
  
  vars_no_interact[ , ,v3[k]] = conv(vars[ , ,v3[k]])
  
  print(paste0('matrix with interactions vars_no_interact[  , ,', v3[k],']'))
  print(vars_no_interact[ , ,v3[k]])
}



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

require('maps')

par(mfrow=c(1,3))
map('usa')
points(x=fia_db$LON[plots_with_tree], y=fia_db$LAT[plots_with_tree])

map('usa')
points(x=fia_db$LON[plots_no_interact], y=fia_db$LAT[plots_no_interact])

map('usa')
points(x=fia_db$LON[plots_interact], y=fia_db$LAT[plots_interact])

