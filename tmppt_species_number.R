#install.packages ('rgdal')
code_path = 'E:/Keefe/Code' # change to where the code is
data_path = 'E:/Keefe/Data' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)

require('raster')
#require('rgdal')

species_number = 212

source('fia_extract.R')
fia_db = fia_extract(p='C:\Users\Olga Rumyantseva\Documents\R Files', species_to_report = species_number)

source('clim_extract.R')
tmpppt = extract_tmpppt (fia_db, p='E:/Keefe/Data/') 

saveRDS(tmpppt, file = paste0('tmpppt_',species_number,'.rds'))
