#OLS estimate of AR1 coefficient
#from Orcutt and Winokur, 1969, Econometrica, 37:1,1-14

OLSAR1 = function(OL)
{
  sumNom = 0
  sumDenom = 0
  Nobs = length(OL)
  ave1 = 0
  ave2 = 0

  for (i in 2:Nobs)
  {
    ave1 = ave1 + OL[i]
    ave2 = ave2 + OL[i - 1]
  }
  ave1 = ave1 / (Nobs - 1)
  ave2 = ave2 / (Nobs - 1)
  for (i in 2:Nobs)
  {

    sumNom = sumNom + (OL[i] - ave1) * (OL[i - 1] - ave2)
    sumDenom = sumDenom + (OL[i - 1] - ave2) * (OL[i - 1] - ave2)
  }
  
  if (sumDenom > 0)
  {
    OLSAR1v = sumNom / sumDenom
  }

  if (sumDenom <= 0)
  {
    OLSAR1v = 0
  }
  return(OLSAR1v)
}