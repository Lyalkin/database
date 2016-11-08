shapleys = function (vars_selected, var)
{
  clim_vars_selected = c(vars_selected, var)
  
  dimensions = paste0('Bio', clim_vars_selected, 's_bin')
  
  fia_db_with_Bins_marker_T = rep(0, n_T)
  
  fia_db_with_Bins_marker_all = rep(0, n_all)
  
  l = length(vars_selected)
  
  for (i in 1:l) 
    {
    
      if (i <= 10) 
        { fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]*res^(i-1)
          fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]*res^(i-1) }
      
       else
         { fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]/res^(i-10)
           fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]/res^(i-10) }
    }
  
  score_without = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
  
  # compute the score with the additional var:
  
  i = i + 1
  
  if (i <= 10) 
    { fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]*res^(i-1)
      fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]*res^(i-1) } 
  else 
     { fia_db_with_Bins_marker_T = fia_db_with_Bins_marker_T + fia_db_with_Bins_T[,dimensions[i]]/res^(i-10)
       fia_db_with_Bins_marker_all = fia_db_with_Bins_marker_all + fia_db_with_Bins[,dimensions[i]]/res^(i-10) }
  
  score_with = sum(fia_db_with_Bins_marker_all %in% unique(na.omit(fia_db_with_Bins_marker_T)))
  shapleys = (score_without - score_with) * ( factorial(l) * factorial(18-l) / factorial(19) ) 
return(shapleys)}

