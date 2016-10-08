# Loads the functions defined in `fia_extract.R`
source('fia_extract.R')


fia_db = fia_extract(p ='C:/Users/Olga Rumyantseva/Documents/R Files/Data/', species_to_report = 212)

saveRDS(fia_db, file = paste0('C:/Users/Olga Rumyantseva/Documents/R Files/Code/', '/fia_db_', 212,'.rds'))

# Look up thoses plot where Giant Sequoias are present:
plots_with_giant_sequoia = which(fia_db$REL_BA_212 > 0)

# Plot them on a map:
require('maps')
map('usa')

points(x=fia_db$LON[plots_with_giant_sequoia], y=fia_db$LAT[plots_with_giant_sequoia])

total_basal_area=sum(fia_db$BASAL.AREA[plots_with_giant_sequoia])





