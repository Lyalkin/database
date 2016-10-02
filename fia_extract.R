################################
### Loading the FIA database ###
################################

## Main function to parse the US FIA inventory:
fia_extract = function(p, species_to_report) 
{
  require(doParallel)
  cl <- #detectCores()#makeCluster(4)
  registerDoParallel(cl)
  states_without_AK_HI=c("AL","AZ","AR","CA","CO","CT","DE","FL","GA","ID","IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")
  # states_without_AK_HI=c("ME")
  CHARAC_FIA = foreach(state = states_without_AK_HI, .combine='rbind', .export=c('get_indices', 'read.csv.onlycols')) %do% {
    cat(paste('\nNow loading FIA data for state',state,'...\n'))
    TREES = read.csv.onlycols(paste0(p,'FIA/',state,"_TREE.CSV"),c(
      "PLT_CN", "INVYR", "SPCD", "DIA", "TPA_UNADJ", 'STATUSCD') )
    COND = read.csv.onlycols(paste0(p,'FIA/',state,"_COND.CSV"), c(
      'PLT_CN', 'STDAGE', 'PHYSCLCD', 'CONDID'))
    PLOT = read.csv.onlycols(paste0(p,'FIA/',state,"_PLOT.CSV"), c(
      'CN', 'LAT', 'LON', 'ECOSUBCD', 'ELEV'))
    
    tmp=get_indices(TREES, rankings, species_to_report)
    tmp = merge(tmp,
                COND,
                by.x='ID', by.y='PLT_CN',
                all.x=T)
    tmp = merge(tmp,
                PLOT,
                by.x='ID', by.y='CN',
                all.x=T)
  }
  #stopCluster(cl)
  CHARAC_FIA$STDAGE[which(CHARAC_FIA$STDAGE<=0)] = NA # these should be NA
  CHARAC_FIA$YEAR[which(CHARAC_FIA$STDAGE==9999)] = NA
  
  # get the soil type:
  type = ifelse(CHARAC_FIA$PHYSCLCD > 10 & CHARAC_FIA$PHYSCLCD < 20, 'xeric', '')
  type = paste0(type, ifelse(CHARAC_FIA$PHYSCLCD > 20 & CHARAC_FIA$PHYSCLCD < 30, 'mesic', ''))
  type = paste0(type, ifelse(CHARAC_FIA$PHYSCLCD > 30 & CHARAC_FIA$PHYSCLCD < 40, 'hydric', ''))
  type[type=='NANANA'] = NA
  CHARAC_FIA$PHYSCLCD = type
  
  # simplify the ecoregion notation:
  CHARAC_FIA$ECOSUBCD = gsub('..$','',CHARAC_FIA$ECOSUBCD)
  return(CHARAC_FIA)
}


# useful function to load only some columns from a big .csv files
read.csv.onlycols = function(path, cols)
{
  XXX=read.csv(path, header=T, nrows=5)
  col_mask = rep('NULL',ncol(XXX))
  col_mask[which(colnames(XXX) %in% cols)] = NA
  return(read.csv(path, header=T, colClasses=col_mask))
}



## This is the main function that calculates the plot indices
get_indices = function(x, rankings, sp)
{
  print(5)
  x.indices = which(x$STATUSCD == 1) # only live trees
  colnames(x)[1] = 'UNIQUEPLOTID'
  y=x[x.indices,colnames(x) %in%  c("UNIQUEPLOTID", "INVYR", "SPCD", "DIA", "TPA_UNADJ")]
  
  results = NULL
  years = sort(unique(y$INVYR))
  cat("Compute stand level indices:\n")
  pb = txtProgressBar(style = 3)
  for(year in years) {
    setTxtProgressBar(pb, (year-years[1]) / (years[length(years)]-years[1]))
    z = y[which(y$INVYR==year),]
    plot.number = unique(z$UNIQUEPLOTID);
    for (i in plot.number) {
      relevant = z[which(z$UNIQUEPLOTID==i),];
      relevant = droplevels(relevant)
      basal.area = (pi*(relevant$DIA/100/2)^2) * relevant$TPA_UNADJ # this is in m²/ha
      stand.basal.area = sum(basal.area)
      
      species_i = which(relevant$SPCD == sp)
      species.basal.area = sum((pi*(relevant$DIA[species_i]/100/2)^2) * relevant$TPA_UNADJ[species_i]) / stand.basal.area
      
      results = rbind(results,c(i, year, stand.basal.area, species.basal.area))
    }
  }
  
  colnames(results)=c("ID","YEAR","BASAL.AREA",paste0('REL_BA_',sp))
  
  return(as.data.frame(results))
}