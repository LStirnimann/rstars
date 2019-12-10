#function for Average Mean
#EstAve = estimated average of the range
#Xdev = deviation from the median
#Xweight = weights

WeightedAverage <- function(W,Sigma_s,Huber)
{
  EstAve = mean(W)
  for (k in 1:2)
  {
    SumofWeights = 0
    SumAve = 0
    for (i in 1:length(W))
    {
      Xdev = ((W[i] - EstAve) / sqrt(Sigma_s))
      #determine the weight of the normalized deviation
      if (is.na(Xdev) || Xdev == '')
      {
        break
      }
      if (Xdev == 0)
      {
        Xweight = 1
      }
      else if (Xdev != 0)
      {
        Xweight = min(1,(Huber / abs(Xdev)))
      }
      
      
      #sum weights and weighed values
      SumofWeights = SumofWeights + Xweight
      SumAve = SumAve + Xdev * Xweight
    }
    SumAve = SumAve / SumofWeights
    SumAve = SumAve * sqrt(Sigma_s) + EstAve
    EstAve = SumAve
  }
  WeightedAverage = EstAve
}