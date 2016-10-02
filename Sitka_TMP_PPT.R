code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('maps')
require('raster')
require('sp')

tmpppt = readRDS(file ='tmpppt_98.rds') 
tmpppt_df = as.data.frame(tmpppt)

#BIO1= Annual Mean Temperature
#BIO12= Annual Precipitation
#summary(tmpppt_df)


plot(tmpppt_df$Bio1s, tmpppt_df$Bio12s, type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(-55, 255), ylim=c(45, 3380), cex=0.5)



Ann_Mean_Temp_cut = cut(tmpppt_df$Bio1s, breaks = 10, labels = FALSE)
Ann_Precip_cut = cut(tmpppt_df$Bio12s, breaks = 10, labels = FALSE)

plot(Ann_Mean_Temp_cut, Ann_Precip_cut, type = "p",
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(0, 11), ylim=c(0, 11), cex=0.5)



























