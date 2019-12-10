
IPN4 = function(est,Nsub)
  #Calculates bias-corrected AR(1) coefficient [est] using
  #interative 4-step procedure with corrections that are
  #reverse proportional to the subsample size - nsub
{
  IPN4v = est + (1 / Nsub)
  for (i in 1:3)
  {
    IPN4v = IPN4v + (abs(IPN4v) / Nsub)
  }
  return(IPN4v)
}