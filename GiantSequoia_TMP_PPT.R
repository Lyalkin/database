code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('maps')
require('raster')
require('sp')
require('rgdal')

species_number = 212
fia_db = readRDS(file ='fia_db_212.rds')


source('clim_extract.R')
tmpppt = extract_tmpppt (fia_db, p='C:/Users/Olga Rumyantseva/Documents/R Files/Data/') 
saveRDS(tmpppt, file = paste0(code_path, '/tmpppt_', 212,'.rds'))


tmpppt_df = as.data.frame(tmpppt)
summary(tmpppt_df)


#BIO1= Annual Mean Temperature
#BIO12= Annual Precipitation


plot(tmpppt_df$Bio1s, tmpppt_df$Bio12s, type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(-55, 255), ylim=c(45, 3380), cex=0.5)



Ann_Mean_Temp_cut = cut(tmpppt_df$Bio1s, breaks = 10, labels = FALSE)
Ann_Precip_cut = cut(tmpppt_df$Bio12s, breaks = 10, labels = FALSE)




plot(Ann_Mean_Temp_cut, Ann_Precip_cut, type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)



tmp_ppt_together = cbind(Ann_Mean_Temp_cut, Ann_Precip_cut)
tmp_ppt_together_df = as.data.frame(tmp_ppt_together)
summary(tmp_ppt_together_df)
u = unique(tmp_ppt_together_df)
u1 = na.omit(u)
u1
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
















  
  
  
  