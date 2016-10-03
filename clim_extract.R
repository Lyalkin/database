####################################################
### Cross-link FIA plots with climatic variables ###
####################################################
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Any results obtained with this code should refer to the following publication:
# Liénard, J. and Harrisson, J. and Strigul, N., "U.S. Forest Response to Projected Climate-Related Stress: a Tolerance Perspective", Global Change Biology (2016), doi: 10.1111/gcb.13291
#####

# green to red color palette with 100 steps:
colm100=rev(colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837'))(100))

## main function:
extract_tmpppt = function(fia_db, p)
{
  require('raster')
  require('rgdal')
  tmps = raster(paste0(p,'WORLDCLIM/bio_1.bil'))
  ppts = raster(paste0(p,'WORLDCLIM/bio_12.bil'))
  pptd = raster(paste0(p,'WORLDCLIM/bio_14.bil'))
  coords = cbind(fia_db$LON,
                 fia_db$LAT)
  coords = SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=NAD83"))
  coords = spTransform(coords, crs(tmps))
  tmps = extract(tmps,coords)
  ppts = extract(ppts,coords)
  pptd = extract(pptd, coords)
  return(list(tmps=tmps, ppts=ppts, pptd = pptd))
}

# shows the tolerance of FIA stands in the climatic space
plot_climatic = function(tmps, ppts, val, title_label, cutoff=NA)
{
  xr = sample(length(val))
  plot(tmps[xr]/10, ppts[xr]/10,
       pch=1,
       col=colm100[cut(val[xr],breaks=seq(0,1,le=101),include.lowest = T)],
       xlim=c(-5,30),bty='n',las=1,
       xlab='Temperature (⁰C)',ylab='Precipitation (mm/month)')
  mtext(title_label,font = 2)
}

# helper function to put a colorbar on a raster plot
fudgeit = function(d, colramp, ...){
  require('fields')
  x = seq(0,1,le=length(d))
  y = seq(0,1,le=length(d))
  image.plot(list(x=x,y=y,z=d), col = colramp, legend.only = T, add=T,...)
}
