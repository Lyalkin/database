conv = function(vars)
{
 vars_no_interact = vars

 for (i in 1:10)
 {
  for (j in 1:10)
  {
    if (vars_no_interact[i,j] == 0) 
    {
      s1 = vars_no_interact[i,1]+vars_no_interact[i,2]+vars_no_interact[i,3]+
        vars_no_interact[i,4]+vars_no_interact[i,5]+vars_no_interact[i,6]+
        vars_no_interact[i,7]+vars_no_interact[i,8]+vars_no_interact[i,9]+
        vars_no_interact[i,10]
      
      s2 = vars_no_interact[1,j]+vars_no_interact[2,j]+vars_no_interact[3,j]+
        vars_no_interact[4,j]+vars_no_interact[5,j]+vars_no_interact[6,j]+
        vars_no_interact[7,j]+vars_no_interact[8,j]+vars_no_interact[9,j]+
        vars_no_interact[10,j]
      
      if (s1>0 & s2>0) {vars_no_interact[i,j] = 1}
    }
  }
 }

 return(vars_no_interact)
}