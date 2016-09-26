# Loads the functions defined in `fia_extract.R`
source('fia_extract.R')

# The following line parses the FIA database and extracts the plots info
# You need to replace the value of `p` to the path on your computer where the FIA directory is located
# Look at the definition of the function `fia_extract` from file `fia_extract.R` to see what it does
fia_db = fia_extract(p='C:\Users\Olga Rumyantseva\Documents\R Files', species_to_report = 212)

# Look up thoses plot where Giant Sequoias are present:
plots_with_giant_sequoia = which(fia_db$REL_BA_212 > 0)

# Plot them on a map:
require('maps')
map('usa')
points(x=fia_db$LON[plots_with_giant_sequoia], y=fia_db$LAT[plots_with_giant_sequoia])

total_basal_area=sum(fia_db$BASAL.AREA[plots_with_giant_sequoia])

#Now plot TMP and PPT for Giant Sequoia


