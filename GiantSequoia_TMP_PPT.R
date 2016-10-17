code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('maps')
require('raster')
require('sp')
require('rgdal')

# 212 is for Giant Sequoia
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



# Now cut (make bins) for all NCOL(tmpppt_df) = 19 climatic variables
source('extract_breaks_bins.R')
breaks = extract_breaks(tmpppt_df)
bins = extract_bins(tmpppt_df, breaks)


#tmppt_df_binned is tmppt_df where each column is binned
tmppt_df_binned = as.data.frame(bins[[1]])
names(tmppt_df_binned) = (paste0(names(tmpppt_df)[1],'_bin',collapse = ''))
for (i in 2:NCOL(tmpppt_df))
{
a = as.data.frame(bins[[i]])
names(a) = (paste0(names(tmpppt_df)[i],'_bin',collapse = ''))
tmppt_df_binned = cbind (tmppt_df_binned, a)
}
head(tmppt_df_binned)
dim(tmppt_df_binned)




# we can plot clim variables
plot(bins[[1]], bins[[12]], type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)

plot(bins[[3]], bins[[5]], type = "p",
     xlab="clim var 3 ", ylab="clim var 15",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)



#model with no interaction

#clim_vars_matrix is a matrix with 2 columns from tmppt_df_binned, na omitted, duplicates omitted

clim_vars_matrix = na.omit(unique(as.data.frame(cbind(tmppt_df_binned[,1], tmppt_df_binned[,12]))))
dim(clim_vars_matrix)
head(clim_vars_matrix)


#vars is a matrix having vars[i,j]=1 if clim_vars_matrix has a row "i j" 
#and                     vars[i,j]=0 otherwise

vars = matrix(0, 10, 10)
for (i in 1:dim(clim_vars_matrix)[1])
{vars[clim_vars_matrix[i,1],clim_vars_matrix[i,2]] = 1}
print(vars)


#vars_no_interact obtained from vars matrix making it 'looking convex'
vars_no_interact = vars
for (i in 1:10)
{
  for (j in 1:10)
  {
    if (vars_no_interact[i,j] == 0) 
    {
      if (vars_no_interact[i,1]+vars_no_interact[i,2]+vars_no_interact[i,3]+
          vars_no_interact[i,4]+vars_no_interact[i,5]+vars_no_interact[i,6]+
          vars_no_interact[i,7]+vars_no_interact[i,8]+vars_no_interact[i,9]+
          vars_no_interact[i,10] > 0) vars_no_interact[i,j] = 1
    
      if (vars_no_interact[1,j]+vars_no_interact[2,j]+vars_no_interact[3,j]+
          vars_no_interact[4,j]+vars_no_interact[5,j]+vars_no_interact[6,j]+
          vars_no_interact[7,j]+vars_no_interact[8,j]+vars_no_interact[9,j]+
          vars_no_interact[10,j] > 0) vars_no_interact[i,j] = 1
    }
  }
}
print(vars_no_interact)



# p is a presence_absence vector for FIA plots if we consider model with no interactions

d = dim(tmppt_df_binned)[1]
p = vector(mode = "logical", length = d)  #d=754851, vector has FALSE components initially

c = as.data.frame(cbind(tmppt_df_binned[,1], tmppt_df_binned[,12]))


for (i in 1:d) {p[i] = vars_no_interact[c[i,1], c[i,2]]}


for (i in 1:d) 
{
  if (is.na(p[i]) == TRUE) {p[i] = 0} 
}




















  
  
  
  