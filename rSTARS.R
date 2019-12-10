# R-package of Sequential T-test Analysis of Regime Shifts (rSTARS)


rstars <- function(data.timeseries = PDO, l.cutoff, pValue = 0.05, Huber = 1, Endfunction = F,
                  preWhitening = F, OLS = F, MPK = F, IP4 = F, SubsampleSize = (l + 1) / 3 ,
                  save.data = T,
                  show.plot = T,
                  FilteredData = T, save.path = (choose.dir()), timeseries = T)
{

  #call the functions
  source("Alpha.R")
  source("EqN.R")
  source("IMPK.R")
  source("IPN4.R")
  source("OLS.R")
  source("WeightedAverage.R")
  source("Stars_Citation.R")
  PDO <- read.table("PDO.txt",header = T, dec = ".")




  ####definition of the parameters####

  TS <- data.timeseries
  l <- l.cutoff
  Plots <- show.plot
  Nsub <- SubsampleSize

  #definition of prewhitening

  if (preWhitening == T)
  {
    if (OLS == F & MPK == F &  IP4 == F)
    {
      stop("preWhitening = T specify OLS, MPK or IP4")
    }

  }

  if (preWhitening == F)
  {
    FilteredData = F
    DT = 0
    if (OLS == T | MPK == T |  IP4 == T)
    {
      stop("preWhitening = F")
    }
  }

  if (preWhitening == TRUE)
  {
    RSI_mat = matrix(0, nrow = length(TS[, 1]) - 1, length(TS[1, ]))
    TabTSpw = matrix(0, nrow = length(TS[, 1]) - 1, length(TS[1, ]))
    TSpw = vector(length = length(TS[, 1]) - 1)
    RMean_mat = matrix(0, nrow = length(TS[, 1]) - 1, length(TS[1, ]))
  }

  if (preWhitening == FALSE)
  {
    RSI_mat = matrix(0, nrow = length(TS[, 1]), length(TS[1, ]))
    RMean_mat = matrix(0, nrow = length(TS[, 1]), length(TS[1, ]))
  }


  #### attaching the data set and removing of red noise ####


  for (TIMESERIESindex in 2:length(TS[1, ]))


  {

    X = ts(TS[, TIMESERIESindex])

    N = length(X)

    if (N < l)
    {
      stop("CutOff cannot be > Time series length")
    }

    #test the subsample size (Nsub) limits
    if (Nsub < 5 &  MPK == TRUE)
    {
      warning("The subsample size is too small. Automatically corrected - minimum value = 5")
      Nsub = 5
    }

    if (Nsub < 3 & (IP4 == TRUE | OLS == TRUE))
    {
      Nsub = 3
      warning("The subsample size is too small. Automatically corrected - minimum value = 3")
    }

    if (Nsub > N)
    {
      Nsub = N
    }

    #-------------------------------------------------------------------------------------------------------
    # Use prewhitening to remove red noise x(t) = x(t) - alpha * x(t-1)

    if (OLS == T | MPK == T | IP4 == T)
    {
      alpha = AlphaEstf(X,N, Nsub,MPK,IP4,OLS)
    }

    if (preWhitening == TRUE)
      #use of prewhitening to remove red noise x(t)=x(t)-alpha*x(t-1)
    {
      for (i in 2:length(X))
      {
        TSpw[i - 1] = X[i] - (alpha * X[(i - 1)])
      }
      X = TSpw
      TabTSpw[, TIMESERIESindex] = TSpw
    }

    #===================#
    ####  STARS 3.2  ####
    #===================#

    #freedom degree
    df = 2 * l - 2

    #two tailed test
    t_stu = abs(qt(pValue / 2, df))

    #Variance and Sigma calcualation for DIFF formula
    A = var(X[1:l])


    for (i in 2:(length(X) - l + 1))
    {
      B = var(X[i:(i + l - 1)])
      A = rbind(A, B)
    }

    #Sigma square
    Sigma_s = mean(A)

    #between mean values of two subsequent regimes that would be statistically
    #significant according to the Student’s t-test
    diff = t_stu * sqrt((2 * Sigma_s) / l)

    #====================#
    #     core steps     #
    #====================#


    vRMean = 0
    RSI = seq(0, 0, length.out = length(X))


    R1 = X[1:l]
    RegimeMean = WeightedAverage(R1,Sigma_s,Huber)
    changepoint = 1
    n1 = 0



    for (intYear in 2:length(X))
    {

     if (is.na(RegimeMean) || RegimeMean == '')
      {
        break
      }


      if (Endfunction == T & intYear == (length(X) - l + 1))
      {
        if (preWhitening == F)
        {
          RSI[(length(X) - l + 1):length(X)] == seq(0, 0, length.out = l)
          break
        }

        if (preWhitening == T)
        {
          RSI[(length(X) - l + 1):(length(X) - 1)] == seq(0, 0, length.out = (l - 1))
          break
        }
      }

      if (X[intYear] > (RegimeMean + diff))
      {
        sumofWeights = 0
        cusumUP = 0
        Xdev = 0
        for (t in intYear:(intYear + l - 1))
        {
          if (t > length(X))
          {
            if (sumofWeights > 0)
            {
              break
            }
          }

          Xdev = (X[t] - RegimeMean - diff) / sqrt(Sigma_s)

          #determine the weight of the normalized deviation
          if (Xdev == 0)
          {
            Xweight = 1
          }

          else if (Xdev != 0)
          {
            Xweight = min(1, (Huber / abs(Xdev)))
          }

          #sum weights and weighed values
          sumofWeights = sumofWeights + Xweight
          cusumUP = cusumUP + (Xdev * Xweight)

          #check if cusum turns zero
          if (cusumUP < 0)
          {
            cusumUP = 0
            break
          }
        }
        cusumUP = cusumUP / sumofWeights

        RSI[intYear] = cusumUP
      }

      else if (X[intYear] < (RegimeMean - diff))
      {
        sumofWeights = 0
        cusumDown = 0
        Xdev = 0
        for (t in intYear:(intYear + l - 1))
        {
          if (t > length(X))
          {
            if (sumofWeights > 0)
            {
              break
            }
          }

          Xdev = (X[t] - RegimeMean + diff) / sqrt(Sigma_s)
          #determine the weight of the normalized deviation
          if (Xdev == 0)
          {
            Xweight = 1
          }
          else if (Xdev != 0)
          {
            Xweight = min(1, (Huber / abs(Xdev)))
          }

          #sum weights and weighed values
          sumofWeights = sumofWeights + Xweight
          cusumDown = cusumDown + (Xdev * Xweight)

          #check if cusum turns zero
          if (cusumDown > 0)
          {
            cusumDown = 0
            break
          }
        }
        cusumDown = cusumDown / sumofWeights
        RSI[intYear] = cusumDown
      }


      else if (RegimeMean - diff <= X[intYear] &
               X[intYear] <= RegimeMean + diff)
      {
        RSI[intYear] = 0
      }

      #check for the situation when the test is not over for the last
      #change point, but we are too close to the end of the time series
      if (abs(RSI[intYear] > 0 & intYear > (length(X) - l + 1)))
      {
        break
      }
      #------------------------------------------------------------------#

      if (RSI[intYear] == 0)
        #intYear is not a new changepoint
      {
        if ((changepoint + l) <= intYear)
        {
          #recalculate regime mean and Diff
          #currently Diff remains constant for the entire process /series
          n1 = intYear - changepoint + 1
          for (n in 1:n1)
          {
            R1[n] = X[changepoint + n - 1]
          }
          RegimeMean = WeightedAverage(R1,Sigma_s,Huber)
        }
      }


      if (RSI[intYear] != 0)
        #regime shift is detected
        #intYear is a new changepoint
      {
        changepoint = intYear
        #recalculate regime mean and Diff
        #currently Diff remains constant for the entire process /series}
        R1 = 0
        for (n in 1:l)
        {
          R1[n] = X[changepoint + n - 1]
        }
        RegimeMean = WeightedAverage(R1,Sigma_s,Huber)

      }
    }

    #Series of RegimeMeans
    if (FilteredData == T)
    {
      S = 1
      for (i in 1:length(RSI))
      {
        if (RSI[i] != 0)
        {
          E = (i - 1)
          MeanRegime = WeightedAverage(X[S:E],Sigma_s,Huber)
          vRMean1 = rep(MeanRegime, length(X[S:E]))
          vRMean = c(vRMean, vRMean1)
          S = i
        }
        if (i == length(RSI))
        {
          E = (length(RSI))
          MeanRegime = WeightedAverage(X[S:E],Sigma_s,Huber)
          vRMean1 = rep(MeanRegime, length(X[S:E]))
          vRMean = c(vRMean, vRMean1)
        }
      }
    }

    if (FilteredData == F)
    {
      X1 = TS[, TIMESERIESindex]
      S = 1
      for (i in 1:length(RSI))
      {
        if (RSI[i] != 0)
        {
          E = (i - 1)
          MeanRegime = WeightedAverage(X1[S:E],Sigma_s,Huber)
          vRMean1 = rep(MeanRegime, length(X1[S:E]))
          vRMean = c(vRMean, vRMean1)
          S = i
        }
        if (i == length(RSI))
        {
          E = (length(RSI))
          MeanRegime = WeightedAverage(X1[S:E],Sigma_s,Huber)
          vRMean1 = rep(MeanRegime, length(X1[S:E]))
          vRMean = c(vRMean, vRMean1)
        }
      }
    }

    vRMean = vRMean[-1]
    RSI_mat[, TIMESERIESindex] = RSI
    RMean_mat[, TIMESERIESindex] = vRMean

  }

  ####Saving tables of regimes avarege (tsMean.txt), RSI.txt and Filtered time series (Filredts.txt)####

  colnames(RMean_mat) <- colnames(data.timeseries)
  colnames(RSI_mat) <- colnames(data.timeseries)

  if (preWhitening == T)
  {
    zeri = seq(0, 0, length.out = length(TS[1, ]))
    RSI_mat = rbind(zeri, RSI_mat)

    empties = rep(NA, length(TS[1, ]))
    RMean_mat = rbind(empties, RMean_mat)


    TabTSpw = rbind(empties, TabTSpw)
    colnames(TabTSpw) <- colnames(data.timeseries)
    TabTSpw_save = TabTSpw

    if (save.data == TRUE)
    {TabTSpw_save[,1] <- data.timeseries[,1]
    path_temp = paste(toString(save.path),"\\Filteredts.txt",sep="")
    write.table(TabTSpw_save, file = path_temp,row.names = FALSE)}
  }

  RMean_mat_save = RMean_mat
  RSI_mat_save = RSI_mat


  if (save.data == TRUE)
  {RMean_mat_save[,1] <- data.timeseries[,1]
  RSI_mat_save[,1] <- data.timeseries[,1]

  path_temp = paste(toString(save.path),"\\tsMean.txt",sep="")
  write.table(RMean_mat_save, file = path_temp,row.names = FALSE)

  path_temp = paste(toString(save.path),"\\RSI.txt",sep="")
  write.table(RSI_mat_save, file = path_temp,row.names = FALSE)
  }





  #### PLOTS ####
  if (Plots == TRUE)
  {
    if (timeseries == TRUE)
    {
    require(xts)
    Time <- as.POSIXct(TS[,1], optional = TRUE, origin = TS[1,1])

    if (preWhitening == F)
    {
      for (i in 2:length(TS[1, ]))
      {
        table_plot = cbind(as.data.frame(RMean_mat)[,i],TS[,i])
        colnames(table_plot) <- c(colnames(as.data.frame(RMean_mat)[i]),colnames(TS[i]))
        p <- plot(
          xts(table_plot,order.by = Time),
          col = c("#EB4C43","#3081B5"),
          #lwd = c(1, 2),
          main = paste("Regime Shift detection in",colnames(TS[i]),"using STARS")
        )
        addLegend('topright',
                  legend.names = c("Time series", "Regimes"),
                  col = c("#3081B5", "#EB4C43"),
                  lwd = c(2, 2)
        )

        print(p)

        rsi <- abs(data.frame(RSI_mat, row.names = NULL))
        nome <- colnames(rsi[i])
        rsi = xts(rsi[,i], order.by = Time)
        colnames(rsi) <- nome
        q <- plot(rsi, type = "h",
             main = paste("Regime shift index values for", colnames(rsi)))
        print(q)

      }
    }
    else if (preWhitening == T)
    {
      if (FilteredData == F)
      {
        for (i in 2:length(TS[1, ]))
        {
          table_plot = cbind(as.data.frame(RMean_mat)[,i],TS[,i])
          colnames(table_plot) <- c(colnames(as.data.frame(RMean_mat)[i]),colnames(TS[i]))
          p <- plot(
            xts(table_plot,order.by = Time),
            col = c("#EB4C43","#3081B5"),
            #lwd = c(1, 2),
            main = paste("Regime Shift detection in",colnames(TS[i]),"using STARS")
          )
          addLegend('topright',
                    legend.names = c("Time series", "Regimes"),
                    col = c("#3081B5", "#EB4C43"),
                    lwd = c(2, 2)
          )

          print(p)

          rsi <- abs(data.frame(RSI_mat, row.names = NULL))
          nome <- colnames(rsi[2])
          rsi = xts(rsi[,2], order.by = Time)
          colnames(rsi) <- nome

          q <- plot(rsi, type = "h",
               main = paste("Regime shift index values for", colnames(rsi)))
          print(q)
        }
      }
      else if (FilteredData == T)
      {


        for (i in 2:length(TS[1, ]))
        {
          table_plot = cbind(as.data.frame(RMean_mat)[,i],as.data.frame(TabTSpw)[,i],TS[,i])
          colnames(table_plot) <- c(colnames(as.data.frame(RMean_mat)[i]),colnames(as.data.frame(TabTSpw)[i]),colnames(TS[i]))
          require(xts)


         p <- plot(
            xts(table_plot,order.by = Time),
            col = c( "#EB4C43","#3081B5","#ECDA61"),
            #lwd = c(1, 1, 2),
            main = paste("Regime Shift detection in",colnames(TS[i]),"using STARS")
          )
          addLegend('topright',
                    legend.names = c("Time series", "Filtered ts", "Regimes"),
                 col = c("#ECDA61", "#3081B5", "#EB4C43"),
                 lwd = c(2, 2, 2)
          )


          print(p)


          rsi <- abs(data.frame(RSI_mat, row.names = NULL))
          nome <- colnames(rsi[i])
          rsi = xts(rsi[,i], order.by = Time)
          colnames(rsi) <- nome

          q <- plot(rsi, type = "h",
               main = paste("Regime shift index values for", colnames(rsi)))

          print(q)

        }
      }
    }
  }

    if (timeseries == FALSE)
    {
      tTS = ts(TS)
      tRMean_mat = ts(RMean_mat)

      if (preWhitening == F)
      {
        for (i in 2:length(TS[1, ]))
        {
          ts.plot(
            tTS[, i],
            tRMean_mat[, i],
            col = c("blue", "red"),
            lwd = c(1, 2),
            xlim = c(0, (length(tTS[, i]) + (
              0.3 * length(tTS[, i])
            )))
          )
          legend(
            x = (length(tTS[, i]) + (0.05 * length(tTS[, i]))),
            y = (max(tTS[, i]) - ((
              max(tTS[, i]) - mean(tTS[, i])
            ) / 2)),
            c("Observed data", "Regimes"),
            col = c("blue", "red"),
            lwd = c(1, 2)
          )
        }
      }
      else if (preWhitening == T)
      {
        if (FilteredData == F)
        {
          for (i in 2:length(TS[1, ]))
          {
            ts.plot(
              tTS[, i],
              tRMean_mat[, i],
              col = c("blue", "red"),
              lwd = c(1, 2),
              xlim = c(0, (length(tTS[, i]) + (
                0.3 * length(tTS[, i])
              )))
            )
            legend(
              x = (length(tTS[, i]) + (0.05 * length(tTS[, i]))),
              y = (max(tTS[, i]) - ((
                max(tTS[, i]) - mean(tTS[, i])
              ) / 2)),
              c("Observed data", "Regimes"),
              col = c("blue", "red"),
              lwd = c(1, 2)
            )
          }
        }
        else if (FilteredData == T)
        {
          for (i in 2:length(TS[1, ]))
          {
            ts.plot(
              tTS[, i],
              TabTSpw[, i],
              RMean_mat[, i],
              col = c("grey", "blue", "red"),
              lwd = c(1, 1, 2),
              xlim = c(0, (length(tTS[, i]) + (
                0.3 * length(tTS[, i])
              )))
            )
            legend(
              x = (length(tTS[, i]) + (0.05 * length(tTS[, i]))),
              y = (max(tTS[, i]) - ((
                max(tTS[, i]) - mean(tTS[, i])
              ) / 2)),
              c("Observed data", "Filtered ts", "Regimes"),
              col = c("grey", "blue", "red"),
              lwd = c(1, 1, 2)
            )
          }
        }
      }
    }

    }

  cat("STARS has completed the analysis. Please, run stars_citation() for references")
}

