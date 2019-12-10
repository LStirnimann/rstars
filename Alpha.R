# ALPHA FUNCTION

# est1 = OLS estimate of aplha

AlphaEstf = function(x,N, Nsub,MPK,IP4,OLS)
{
  source("IMPK.R")
  source("IPN4.R")
  source("OLS.R")
  
  #Subsampling and OLS process

  ss = vector(length = (N - Nsub + 1))
  yy = vector(length = Nsub)
  
  for (i in  1:(N - Nsub + 1))
  {
    for (k in 1:Nsub)
    {
      yy[k] = x[i + k - 1]
    }
    ss[i] = OLSAR1(yy)
  }

  est1 = median(ss)
  
  #bias correction of OLS if requested
  if (MPK == TRUE)
  {
    AlphaEst = IMPK(est1,Nsub)
  }
  else if (IP4 == TRUE)
  {
    AlphaEst = IPN4(est1,Nsub)
  }
  else if (OLS == TRUE)
  {
    AlphaEst = est1
  }
  
  
  #limits
  if (AlphaEst < 0)
  {
    AlphaEst = 0
  }
  if (AlphaEst > 0.99)
  {
    AlphaEst = 0.99
  }
  return (AlphaEst)
}
