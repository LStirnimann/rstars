EqN = function (Ns,alpha)
{
  summa = 0
  #calculates the equivalent sample size
  #as in von Storch and Zwiers (1999, p. 115)
  #Note: The formula for the equiv. sample size published by
  #Zwiers and Storch (1995)in J. Climate is somewhat different
  #Ns original sample size
  #alpha - AR1 coefficient defined for the entire module
  
  for (i in 1:(Ns - 1))
  {
    summa = summa + (1 - i / Ns) * alpha ^ i
  }
  EqNv = Ns / (1 + summa)
  
  #just in case
  if (EqNv <= 2)
  {
    EqNv = 2
  }
  if (EqNv > Ns)
  {
    EqNv = Ns
  }
  return(EqNv)
}