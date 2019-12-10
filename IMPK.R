IMPK = function(est,Nsub)
{
  #Calculates bias-corrected AR(1) coefficient using the formula of
  #Marriott-Pope and Kendall(see Orcutt and Winokur, 1969, Econometrica, 37:1,1-14)
  #est is the OLS estimate of AR(1)
  #Nsub - sample size declared globaly

  if (Nsub > 4)
  {
    IMPKv = ((Nsub - 1) * est + 1) / (Nsub - 4)
  }
  else
    #should not be here
  {
    IMPKv = est
  }
  return(IMPKv)
}