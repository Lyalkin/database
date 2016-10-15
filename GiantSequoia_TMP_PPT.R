code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('maps')
require('raster')
require('sp')
require('rgdal')

# 212 is for Giant Sequoias
species_number = 212

# extraction of all the plots in the FIA DB and computation of the relative BA of the species
# source('fia_extract.R')
# fia_db = fia_extract(p=data_path, species_number)
#saveRDS(fia_db, file = paste0(code_path, '/fia_db_', species_number, '.rds'))
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


par(mfrow=c(1,2))

plot(tmpppt_df$Bio1s[plots_with_species], tmpppt_df$Bio12s[plots_with_species], type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(-55, 255), ylim=c(45, 3380), cex=0.5)

breaks_tmp = seq(from = min(tmpppt_df$Bio1s, na.rm=T),
                 to = max(tmpppt_df$Bio1s, na.rm=T),
                 le=11)

Ann_Mean_Temp_cut = cut(tmpppt_df$Bio1s, breaks = breaks_tmp, labels = FALSE)


breaks_ppt = seq(from = min(tmpppt_df$Bio12s, na.rm=T),
                 to = max(tmpppt_df$Bio12s, na.rm=T),
                 le=11)
Ann_Precip_cut = cut(tmpppt_df$Bio12s, breaks = breaks_ppt, labels = FALSE)


plot(Ann_Mean_Temp_cut[plots_with_species], Ann_Precip_cut[plots_with_species], type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)

par(mfrow=c(1,1))

## Now cut (make bins) for all NCOL(tmpppt_df) = 19 climatic variables
source('extract_breaks_bins.R')
breaks = extract_breaks(tmpppt_df)
bins = extract_bins(tmpppt_df)


#tmppt_df_binned is tmppt_df where each column is binned
tmppt_df_binned = as.data.frame(bins[[1]])
names(tmppt_df_binned) = (paste0(names(tmpppt_df)[1],'_bin',collapse = ''))
for (i in 2:19)
{
a = as.data.frame(bins[[i]])
names(a) = (paste0(names(tmpppt_df)[i],'_bin',collapse = ''))
tmppt_df_binned = cbind (tmppt_df_binned, a)
}
head(tmppt_df_binned)
dim(tmppt_df_binned)





plot(bins[[1]], bins[[12]], type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)

plot(bins[[3]], bins[[5]], type = "p",
     xlab="clim var 3 ", ylab="clim var 15",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)



#model with no interaction----------------------------------------------------------

#clim_vars_matrix is a matrix with 2 columns - clim variables bins
var_cbind = cbind(bins[[1]], bins[[12]])
var_cbind_df = as.data.frame(var_cbind)
var_cbind_df_unique = unique(var_cbind_df)
clim_vars_matrix = na.omit(var_cbind_df_unique)
dim(clim_vars_matrix)
head(clim_vars_matrix)


#vars is a matrix having vars[i,j]=1 if clim_vars_matrix has a row "i j" 
#and                     vars[i,j]=0 otherwise

vars = matrix(0, 10, 10)
for (i in 1:dim(clim_vars_matrix)[1])
{vars[clim_vars_matrix[i,1],clim_vars_matrix[i,2]] = 1}
print(vars)

#vars_no_interact obtained from vars matrix
for (i in 1:10)
{
  for (j in 1:10)
  {
    if (vars[i,j] == 0) 
    {
      if (vars[i,1]+vars[i,2]+vars[i,3]+vars[i,4]+vars[i,5]+
        vars[i,6]+vars[i,7]+vars[i,8]+vars[i,9]+vars[i,10] > 0) vars[i,j] = 1
    
      if (vars[1,j]+vars[2,j]+vars[3,j]+vars[4,j]+vars[5,j]+
        vars[6,j]+vars[7,j]+vars[8,j]+vars[9,j]+vars[10,j] > 0) vars[i,j] = 1
    }
  }
}
print(vars)



















  
  
  
  