

BinGeneration = function(tmppptDF,res=10){
  slist = list()
  sdata.frame = data.frame()
  res = as.numeric(res)
  lengthTEMP = NCOL(tmppptDF)
  for (i in 1:NCOL(tmppptDF)) {
    namTEMP = (paste0(names(tmppptDF)[i],'_bin',collapse = ''))
    breaksTEMP = seq(from=min((tmppptDF[[i]]), na.rm=T), to=max((tmppptDF[[i]]), na.rm=T), length.out = res+1)
    value = cut(tmppptDF[[i]], breaks=breaksTEMP, labels=F)
    slist[[namTEMP]] = value
  }
  sdata.frame = as.data.frame(slist)
  return(sdata.frame)
}


BinGeneration_Raster = function(tmppptDF,res=10,extent,extents_name,clim_vars_in_use){
  load('mask_usa.RData')
  slist = list()
  sdata.frame = data.frame()
  res = as.numeric(res)
  count = 1
  for (i in 1:length(clim_vars_in_use)) {
    if (clim_vars_in_use[i] == T){
      namTEMP = (paste0(names(tmppptDF)[count],'_',extents_name,'_bin',collapse = ''))
      breaksTEMP = seq(from=min((tmppptDF[[count]]), na.rm=T), to=max((tmppptDF[[count]]), na.rm=T), length.out = res+1)
      count = count + 1
      rastervalue = crop(raster(paste0(data_path,'WORLDCLIM/bio_',i,'.bil')), extent)
      values(rastervalue)[which(is.na(values(mask_usa)))] = NA
      # plot(rastervalue)
      # title(names(clim_vars_in_use[i]))
      slist[[namTEMP]] = cut(values(rastervalue), breaks=breaksTEMP, labels=FALSE)
    }
  }
  #sdata.frame = as.data.frame(slist)
  return(slist)
  #return(sdata.frame)
}

TrimRaster = function(tmppptDF,res=10,extent,extents_name,clim_vars_in_use,newraster){
  slist = list()
  sdata.frame = data.frame()
  res = as.numeric(res)
  count = 1
  require('raster')
  newextent <- calc(newraster, fun=function(x){replace(x, x==0,NA)})
  for (i in 1:length(clim_vars_in_use)) {
    if (clim_vars_in_use[i] == T){
      namTEMP = (paste0(names(tmppptDF)[count],'_',extents_name,'_bin',collapse = ''))
      breaksTEMP = seq(from=min((tmppptDF[[count]]), na.rm=T), to=max((tmppptDF[[count]]), na.rm=T), length.out = res+1)
      count = count + 1
      rastervalue = crop(raster(paste0(data_path,'WORLDCLIM/bio_',i,'.bil')), extent)
      values(rastervalue)[which(is.na(values(newextent)))] = NA
      # plot(rastervalue)
      # title(names(clim_vars_in_use[i]))
      slist[[namTEMP]] = cut(values(rastervalue), breaks=breaksTEMP, labels=FALSE)
    }
  }
  #sdata.frame = as.data.frame(slist)
  return(slist)
  #return(sdata.frame)
}

#Generating a New tmpppt

trim_tmpppt = function(fia_db,p,clim_vars_in_use,newraster){
  require('raster')
  require('rgdal')
  slist = list()
  firstTrue = T
  newextent = calc(newraster, fun=function(x){replace(x, x==0,NA)})
  for (i in 1:length(clim_vars_in_use)) {#Creates the raster of each climate variable we are using e.g. BIO1s 
    if (clim_vars_in_use[i] == T) {
      namTEMP = (paste0("Bio",i,'s',collapse = ''))
      value = crop(raster(paste0(p,'WORLDCLIM/bio_',i,'.bil')), extent_limit)
      values(value)[which(is.na(values(newextent)))] = NA
      if (firstTrue == T) {
        coords = cbind(fia_db$LON,
                       fia_db$LAT)
        coords = SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=NAD83"))
        coords = spTransform(coords, crs(value))
        firstTrue = F
      }
      value = extract(value,coords)
      slist[[namTEMP]] = value
    }
  }
  #tmps = raster(paste0(p,'WORLDCLIM/bio_1.bil')) #These three replaced by the bio,i,stemp version of each clim var
  #ppts = raster(paste0(p,'WORLDCLIM/bio_12.bil')) 
  #pptd = raster(paste0(p,'WORLDCLIM/bio_14.bil'))
  #tmps = extract(tmps,coords)
  #ppts = extract(ppts,coords)
  #pptd = extract(pptd, coords)
  return(slist)
}









###Generating the raster for each of the best predictions



#3D#######

bestpredictionraster3d = function(possible_combinations_Nd,number_of_predictions_Nd,tmpppt_Bins,tmpppt_raster_bins,res){
  combinations_Nd_potential=matrix(nrow = number_of_predictions_Nd, ncol = (ncol(possible_combinations_Nd)+1))
  predictions_matrix = raster()
  for (i in 1:number_of_predictions_Nd) {
    source('nd_models.R')
    modelND = generate_model_nD(fia_db_with_Bins, res, data_name='PRESENCE_261', 
                                dim1 = paste0('Bio',possible_combinations_Nd[i,2],'s_bin'), 
                                dim2 = paste0('Bio',possible_combinations_Nd[i,3],'s_bin'),
                                dim3 = paste0('Bio',possible_combinations_Nd[i,4],'s_bin'))
    predictionNd = readRDS(file='prediction2d_261.rds')
    values(predictionNd) = NA
    names(predictionNd) = 'PredictedND'
    values(predictionNd) = extract_model_nD(modelND, tmpppt_raster_bins[[possible_combinations_Nd[i,2]]],
                                            tmpppt_raster_bins[[possible_combinations_Nd[i,3]]],
                                            tmpppt_raster_bins[[possible_combinations_Nd[i,4]]])
    plot(predictionNd, main= paste0('Potential Distribution of Eastern Hemlock_Bio',
                                    possible_combinations_Nd[i,2],
                                    '&Bio',possible_combinations_Nd[i,3],
                                    '&Bio',possible_combinations_Nd[i,4]))
    map('state', add=T)
    #points(fia_db$LON, fia_db$LAT, cex=fia_db$PRESENCE_261)
    source('niche_fcts.R')
    #computed_potential = compute_potential(predictionNd)
    #combinations_Nd_potential[i,] = c(computed_potential,
                                      # possible_combinations_Nd[i,1],
                                      # possible_combinations_Nd[i,2],
                                      # possible_combinations_Nd[i,3])
    #combinations_Nd_potential = as.data.frame(combinations_Nd_potential)
  }
  # combinations_Nd_potential = t(as.matrix.default(combinations_Nd_potential))
  # combinations_Nd_potential = combinations_Nd_potential[order(combinations_Nd_potential[,1]),]
  return(predictionNd)
}


