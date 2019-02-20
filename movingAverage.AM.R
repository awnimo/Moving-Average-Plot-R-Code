movingAverage.AM <- function( df, winSize=200, step=40){
  ################################################################################################
  # This function will run sliding window of a moving average method given a window size and steps
  #     input:
  #     ------
  #         df - a two column data frame of a response and terms, where response is the (numeric)
  #              response vector and terms is a series of terms which specifies a linear predictor
  #              for response. Example, response = log2FC, predictor = geneLength
  #         winSize - the desired window size; Integer > 1; default 200
  #         step - the desired distance between windows; Integer <= winSize; default 40

  #     output:
  #     -------
  #         a LIST of 2 elements:
  #           I - a 4 column data frame:
  #             predMean - the mean of the predictor terms within a window
  #             respAvg - the mean of the response within a window (the moving average)
  #             ucl, lcl - the 95% CI upper and lower levels
  #           II - a plotting object

  #     value:
  #     ------
  #         the output can be assigned to a variable to be used in the environment by calling
  #         the first element in the list[[1]] to get the data frame with the values, and by
  #         calling the second element in the list[[2]] to initiate the plot
  ################################################################################################

  colnames(df) <- c("response","predictor")
  meanValue <- NULL
  l <- seq(1,length(df$response),step)
  for(i in l){
    print(i)
    if(i+winSize-1 > length(df$response)){
      i1 = length(df$response) }
    else i1 = i+winSize-1
    # mean
    respAvg <- mean(df$response[i:i1])
    # 95% CI
    ucl <- respAvg + 1.96 * sd(df$response[i:i1]/sqrt(length(df$response[i:i1])))
    lcl <- respAvg - 1.96 * sd(df$response[i:i1]/sqrt(length(df$response[i:i1])))
    # average gene length in bin
    predMean <- mean(df$predictor[i:i1])
    meanValue.i <- cbind(predMean,respAvg,ucl,lcl)
    meanValue <- rbind(meanValue, meanValue.i)
  }
  meanValue <- as.data.frame(meanValue)
  # make plot
  library(ggplot2)
  p = ggplot(df, aes(y=response, x=predictor))
  P <- p + scale_x_log10(breaks=10^(seq(-1,10)),labels=10^(seq(-1,10))) +
    geom_point(color="grey64", alpha=1/2) +
    geom_abline(intercept = 0, slope = 0, lwd = 1, col = "#009E73", alpha=1/2) +
    annotation_logticks(sides="b") +
    geom_smooth(aes(x=predMean, y=respAvg, ymin=lcl, ymax=ucl), data=meanValue, stat="identity", color="red", fill="red") +
    geom_rug(position="jitter", size=0.05, color="magenta", alpha=1/8, sides="tr") +
    xlab("predictor") + ylab("response") +
    ggtitle(bquote( atop("Moving Average", paste("Window Size ",.(winSize), " Step ", .(step) ))))
  return(list(meanValue, P))
}
