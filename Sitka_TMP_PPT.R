code_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Code' # change to where the code is
data_path = 'C:/Users/Olga Rumyantseva/Documents/R Files/Data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('maps')
require('raster')
require('sp')


BIO1 = T ; names(BIO1) = "Annual Mean Temperature"
BIO2 = T ; names(BIO2) = "Mean Diurnal Range (Mean of monthly (max temp - min temp))"
BIO3 = T ; names(BIO3) = "Isothermality (BIO2/BIO7) (* 100)"
BIO4 = T ; names(BIO4) = "Temperature Seasonality (standard deviation *100)"
BIO5 = T ; names(BIO5) = "Max Temperature of Warmest Month"
BIO6 = T ; names(BIO6) = "Min Temperature of Coldest Month"
BIO7 = T ; names(BIO7) = "Temperature Annual Range (BIO5-BIO6)"
BIO8 = T ; names(BIO8) = "Mean Temperature of Wettest Quarter"
BIO9 = T ; names(BIO9) = "Mean Temperature of Driest Quarter"
BIO10 = T ; names(BIO10) = "Mean Temperature of Warmest Quarter"
BIO11 = T ; names(BIO11) = "Mean Temperature of Coldest Quarter"
BIO12 = T ; names(BIO12) = "Annual Precipitation"
BIO13 = T ; names(BIO13) = "Precipitation of Wettest Month"
BIO14 = T ; names(BIO14) = "Precipitation of Driest Month"
BIO15 = T ; names(BIO15) = "Precipitation Seasonality (Coefficient of Variation)"
BIO16 = T ; names(BIO16) = "Precipitation of Wettest Quarter"
BIO17 = T ; names(BIO17) = "Precipitation of Driest Quarter"
BIO18 = T ; names(BIO18)= "Precipitation of Warmest Quarter"
BIO19 = T ; names(BIO19) = "Precipitation of Coldest Quarter"

clim_vars_in_use = c(BIO1,BIO2,BIO3,BIO4,BIO5,BIO6,BIO7,BIO8,BIO9,
                     BIO10,BIO11,BIO12,BIO13,BIO14,BIO15,BIO16,BIO17,BIO18,BIO19)

tmpppt = readRDS(file='tmpppt_98.rds') 
tmpppt_df = as.data.frame(tmpppt)

#BIO1= Annual Mean Temperature
#BIO12= Annual Precipitation
#summary(tmpppt_df)

plot(tmpppt_df$BIO1, tmpppt_df$BIO12, 
     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
     xlim=c(100,250), ylim=c(1000,1500), cex=0.5)

##tmpppt_df_graph=na.omit(tmpppt_df)
#plot(tmpppt_df_graph$BIO1, tmpppt_df_graph$BIO12, 
#     xlab="Annual Mean Temperature ", ylab="Annual Precipitation",
#     xlim=c(-50,250), ylim=c(50,3380), cex=0.5)


test=which(tmpppt_df$BIOs1 == 250)
dim(test)

















