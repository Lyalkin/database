#makes a list of 19 vectors-breaks for each of 19 clim. variables

extract_breaks = function(tmpppt_df)
{
  breaks = list()
  
  for (j in 1:NCOL(tmpppt_df)){
    
    
    breaks[[j]] = seq(from = min(tmpppt_df[,j], na.rm=T),
                      to = max(tmpppt_df[,j], na.rm=T),
                      le = 11)
  }
  return(breaks)
}

#makes a list of 19 vectors-bins(clim. variables cuts) for each of 19 clim. variables
extract_bins = function(tmpppt_df)
{
  bins = list()
  
  for (j in 1:NCOL(tmpppt_df)){
    
    
    bins[[j]] = cut(tmpppt_df[,j], breaks = breaks[[j]], labels = FALSE)
  }
  return(bins)
}







