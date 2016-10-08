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
# fia_db = fia_extract(p=data_path, species_number) ### RUN THIS CODE TO MAKE SURE THAT IT WORKS!! :)
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

## Todo: write a function that cuts the climatic data
breaks_and_bins = extract_breaks_bins(tmpppt_df)
breaks_bio = breaks_and_bins$breaks_bio
bins_bio = breaks_and_bins$bins_bio

plot(bins_bio[[1]], bins_bio[[12]], type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)

plot(bins_bio[[3]], bins_bio[[15]], type = "p",
     xlab="clim var 3 ", ylab="clim var 15",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)



tmp_ppt_together = cbind(Ann_Mean_Temp_cut_species, Ann_Precip_cut_species)
tmp_ppt_together_df = as.data.frame(tmp_ppt_together)
summary(tmp_ppt_together_df)
u = unique(tmp_ppt_together_df)
u1 = na.omit(u)
print(u1)
dim(u1)

tmp_ppt_no_interact = matrix(0, 10, 10)

for (i in 1:67) 
{tmp_ppt_no_interact[u1[i,1], u1[i,2]] = 1}
print(tmp_ppt_no_interact)

a = tmp_ppt_no_interact
print(a)

for (i in 1:9)
{
  for (j in 10:2) 
  {
    if (a[i,j-1] == 0 & a[i,j] == 1 & a[i+1,j-1] == 1 & a[i+1,j] == 1)  a[i,j-1] = 1
    if (a[i,j-1] == 1 & a[i,j] == 1 & a[i+1,j-1] == 0 & a[i+1,j] == 1)  a[i,j-1] = 1
    } 
}
print(a)


for (i in 1:9)
{
  for (j in 1:9) 
  {
    if (a[i,j] == 1 & a[i,j+1] == 1 & a[i+1,j] == 1 & a[i+1,j+1] == 0)  a[i+1,j+1] = 1
    if (a[i,j] == 1 & a[i,j+1] == 0 & a[i+1,j] == 1 & a[i+1,j+1] == 1)  a[i,j+1] = 1
  }
}
print(a)
















  
  
  
  