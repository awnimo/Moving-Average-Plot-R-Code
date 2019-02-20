movingAvrg.lm <- function(df, win=40, step=200, xAxis.indx=seq(1,dim(df)[1])){
  # df-         input data frame with 2 columns: 1st column: variable, 2nd column: values
  # win-        windows length for averaging
  # step-       windows slides by step
  # xAxis.indx- x-axis average size value for win over a sliding step
  meanValue <- NULL
  l <- seq(1,dim(df)[1],win)
  for(i in l){
    print(i)
    if(i+step-1 > dim(df)[1]){
      i1 = dim(df)[1] }
    else i1 = i+step-1
    # mean
    avg <- mean(df[i:i1,2])
    # 95% CI
    ucl <- avg + 1.96 * sd(df[i:i1,2]/sqrt(length(df[i:i1,2])))
    lcl <- avg - 1.96 * sd(df[i:i1,2]/sqrt(length(df[i:i1,2])))
    # average gene length in bin
    xAxis <- mean(xAxis.indx[i:i1])
    meanValue.i <- cbind(xAxis,avg,ucl,lcl)
    meanValue <- rbind(meanValue, meanValue.i)
  }
  return( meanValue )
}