## isn(): implement Interval-weighted SegNeigh (ISN) Changepoint Detection
## INPUT LIST:
## "data": numerical vector, a vector contains the time series to be detected
## "func": string, the fitting method, now only "mean" and "line" are implemented
## "m": non-negative integer, the number of maximum potential changepoints to be detected
## "pnum": non-negative, the hyper-parameter to control penalty on number of changepoints
## "ploc": non-negative, the hyper-parameter to control penalty on length of intervals for segments
## "efficiency": TRUE or FALSE, the method of approximate ISN to be applied or not
## OUTPUT LIST:
## "location": optimal changepoint location sets for 1:m changepoints condition respectively
## "cost": corresponding cost value for 1:m changepoints condition respectively
## "la": record of location during iteration
## "ca": record of cost during iteration
## "n": length of data
## "bestnum": optimal changepoints' number set with lowest cost
## "bestloc": optimal changepoints' location set with lowest cost
## ¡°approx¡±: TRUE or FALSE, whether the approximate algorithm is applied or not

isn <- function(data, func="mean", m=5, pnum=1, ploc=1, efficiency=FALSE) # Interval-weighted SegNeigh Changepoint Detection
{
  data <- as.numeric(data)
  if(is.vector(data)&&is.atomic(data)&&all(!is.na(data))==FALSE)
    stop("data should be vector.")
  if(m%%1!=0)
    stop("m should be integer.")
  if(m<0)
    stop("m should be non-negative.")
  if(pnum<0)
    stop("pnum should be non-negative.")
  if(ploc<0)
    stop("ploc should be non-negative.")
  if(ploc==0)
  {
    efficiency <- TRUE
    # warning("ploc=0, approximate algorithm is applied (with the same result of ISN)")
  }
  n <- length(data)
  location <- matrix(nrow=m, ncol=m) # location record for m
  cost <- rep(0,m+1) # cost record for m
  cm <- matrix(nrow=n, ncol=n) # initial cost matrix
  # STEP 1
  if(func=="mean") # y=b
  {
    dim <- 1
    for(i in 1:n)
      for(j in i:n)
        cm[i,j] <- sum((data[i:j]-mean(data[i:j]))^2)
  }else if(func=="line"){ # y=ax+b
    dim <- 2
    for(i in 1:n)
      for(j in i:n)
        cm[i,j] <- deviance(lm(data[i:j]~1+c(i:j)))
  }
  cost[1] <- cm[1,n]+pnum*dim*log(n)
  upperbound <- cost[1]*log(n)
  # STEP 2
  if(efficiency==FALSE)
  {
    la <- array(dim=c(m,m,n-1)) # location array
    ca <- array(dim=c(m,m,n-1)) # cost array
  }else if(efficiency==TRUE){
    la2 <- array(dim=c(m+1,n-1)) # location array 2
    ca2 <- array(dim=c(m+1,n-1)) # cost array 2
    cmk <- matrix(nrow=n,ncol=n) # cost matrix for k
    cmk <- cm+ploc*cploc2(n=n,m=m)
    ca2[1,] <- cmk[1,1:(n-1)]
  }
  for(k in 1:m) # dynamic programming
  {
    if(efficiency==FALSE)
    {
      cmk <- matrix(nrow=n,ncol=n) # cost matrix for k
      cmk <- cm+ploc*cploc(n=n,k=k)
      if(k==1)
      {
        cpmin <- cpargmin(n=n,j=n,ca=cmk[1,],cmk=cmk[,n],upperbound=upperbound,pnum=pnum,dim=dim)
        location[k,k] <- cpmin[1]
        cost[k+1] <- cpmin[2]+pnum*(dim*(k+1)+k)*log(n)
      }else if(k>1){
        ca[k,1,] <- cmk[1,1:(n-1)]
        for(kk in 1:(k-1))
          for(j in (1+kk):(n-k+kk))
          {
            cpmin <- cpargmin(n=n,j=j,ca=ca[k,kk,],cmk=cmk[,j],upperbound=upperbound,pnum=pnum,dim=dim)
            la[k,kk+1,j] <- cpmin[1]
            ca[k,kk+1,j] <- cpmin[2]
          }
        cpmin <- cpargmin(n=n,j=n,ca=ca[k,k,],cmk=cmk[,n],upperbound=upperbound,pnum=pnum,dim=dim)
        location[k,k] <- cpmin[1]
        cost[k+1] <- cpmin[2]+pnum*(dim*(k+1)+k)*log(n)
        for(kk in (k-1):1)
          location[k,kk] <- la[k,kk+1,location[k,kk+1]]
      }
    }else if(efficiency==TRUE){
      if(k<m)
        for(j in (k+1):(n-1))
        {
          cpmin <- cpargmin(n=n,j=j,ca=ca2[k,],cmk=cmk[,j],upperbound=upperbound,pnum=pnum,dim=dim)
          la2[k+1,j] <- cpmin[1]
          ca2[k+1,j] <- cpmin[2]
        }
      cpmin <- cpargmin(n=n,j=n,ca=ca2[k,],cmk=cmk[,n],upperbound=upperbound,pnum=pnum,dim=dim)
      location[k,k] <- cpmin[1]
      cost[k+1] <- cpmin[2]+pnum*(dim*(k+1)+k)*log(n)
      if(k>1)
        for(kk in (k-1):1)
          location[k,kk] <- la2[kk+1,location[k,kk+1]]
    }
  }
  # OUTPUT
  bestnum <- which.min(cost)-1
  bestloc <- location[bestnum,]
  bestloc <- subset(bestloc,bestloc!="NA")
  if(efficiency==FALSE)
    isn <- list(location=location,cost=cost,data=data,func=func,la=la,ca=ca,n=n,m=m,pnum=pnum,ploc=ploc,bestnum=bestnum,bestloc=bestloc,approx=efficiency)
  if(efficiency==TRUE)
    isn <- list(location=location,cost=cost,data=data,func=func,la=la2,ca=ca2,n=n,m=m,pnum=pnum,ploc=ploc,bestnum=bestnum,bestloc=bestloc,approx=efficiency)
  class(isn) <- "isn"
  return(isn)
}

## cploc(): compute penalty of interval, no user interaction

cploc <- function(n, k) # compute penalty of interval
{
  cl <- matrix(nrow=n, ncol=n)
  for(i in 1:n)
    for(j in i:n)
      cl[i,j] <- ((k+1)*(j-i+1)/n-1)^2*log(n)
  return(cl)
}

## cploc(): compute penalty of interval for approximate ISN, no user interaction

cploc2 <- function(n, m, plot=FALSE) # compute penalty of interval
{
  if(plot==FALSE)
    cl <- matrix(0,nrow=n,ncol=n)
  else if(plot==TRUE)
    cl <- matrix(nrow=n,ncol=n)
  intv <- floor(n/(m+1))
  for(i in 1:n)
    for(j in i:min(i+intv-1,n))
      cl[i,j] <- ((m+1)*(j-i+1)/n-1)^2*log(n)
  return(cl)
}

## cpargmin(): compute argmin for v, no user interaction

cpargmin <- function(n, j, ca, cmk, upperbound, pnum, dim) # compute argmin for v
{
  tmp <- rep(upperbound,j-1)
  for(v in 1:(j-1))
    tmp[v] <- ca[v]+cmk[v+1]
  tmpmin <- which.min(tmp)
  cpmin <- c(tmpmin,tmp[tmpmin])
  return(cpmin)
}

## print.isn(): output summary

print.isn <- function(x, ...)
{
  print("Best changepoints' number:")
  print(x$bestnum)
  print("Best changepoints' locations:")
  print(x$bestloc)
  print("Changepoints location records:")
  print(x$location)
  print("Changepoints costs:")
  print(x$cost)
  print("Fitting method:")
  print(x$func)
  hp <- matrix(c(x$m,x$pnum,x$ploc),nrow=1,ncol=3)
  colnames(hp) <- c("upper-number","penalty-number","penalty-location")
  print(hp)
  if(x$approx==TRUE)
    print("Approximate algorithm is applied")
}

## plot.isn(): output figures
## INPUT LIST:
## "output": numerical value acceptable for 0, 1, 2, 3, 4:
##           0(default): output all 4 types of figures on one plot
##           1: output the input data with fitting model
##           2: output the cost function varied by number of changepoints
##           3: output the squared residual
##           4: output the residual sum of squares
## "xlabel": string, label of x-axis for plot 1
## "ylabel": string, label of y-axis for plot 1

plot.isn <- function(x, output=0, xlabel="t", ylabel="data", ...)
{
  if(which.min(x$cost)==1)
    stop("no change point, alternative hypothesis is rejected.")
  # Computation
  sgnum <- which.min(x$cost)
  cpnum <- sgnum-1
  cploc <- x$location[cpnum,]
  cploc <- subset(cploc,cploc!="NA")
  sgloc <- c(0,cploc,x$n)
  sgmean <- rep(0,sgnum)
  sgr <- rep(0,x$n)
  cf <- matrix(nrow=sgnum,ncol=2)
  for(i in 1:sgnum)
  {
    t1 <- sgloc[i]+1
    t2 <- sgloc[i+1]
    if(x$func=="mean") # y=b
    {
      sgmean[i] <- mean(c(x$data[t1:t2]))
      for(j in t1:t2)
        sgr[j] <- x$data[j]-sgmean[i]
    }else if(x$func=="line"){ # y=ax+b
      lmm=lm(x$data[t1:t2]~1+c(t1:t2))
      cf[i,1] <- lmm$coefficients[[1]]
      cf[i,2] <- lmm$coefficients[[2]]
      for(j in t1:t2)
        sgr[j] <- (x$data[j]-(cf[i,2]*j+cf[i,1]))
    }
  }
  cmrs <- cumsum(sgr^2) # cumulative residual square
  cpcol <- rep("black",x$m+1)
  cpcol[cpnum+1] <- "red"
  # Plot
  if(output==0)
  {
    layout(matrix(c(1,3,2,4),2))
    plot(x=c(1:x$n),y=x$data,xlab=xlabel,ylab=ylabel)
    for(i in 1:sgnum)
    {
      t1 <- sgloc[i]+1
      t2 <- sgloc[i+1]
      if(x$func=="mean") # y=b
        segments(t1,sgmean[i],t2,sgmean[i],col='red')
      else if(x$func=="line") # y=ax+b
        segments(t1,cf[i,2]*t1+cf[i,1],t2,cf[i,2]*t2+cf[i,1],col='red')
    }
    plot(x=c(0:x$m),y=x$cost,type="b",xlab="number of changepoints",ylab="cost",col=cpcol)
    plot(x=c(1:x$n),y=sgr,type="l",xlab="t",ylab="Residual")
    plot(x=c(1:x$n),y=cmrs,type="l",xlab="t",ylab="Cumulative Residual Square")
  }else if(output==1){
    plot(x=c(1:x$n),y=x$data,xlab=xlabel,ylab=ylabel)
    for(i in 1:sgnum)
    {
      t1 <- sgloc[i]+1
      t2 <- sgloc[i+1]
      if(x$func=="mean") # y=b
        segments(t1,sgmean[i],t2,sgmean[i],col='red')
      else if(x$func=="line") # y=ax+b
        segments(t1,cf[i,2]*t1+cf[i,1],t2,cf[i,2]*t2+cf[i,1],col='red')
    }
  }else if(output==2){
    plot(x=c(0:x$m),y=x$cost,type="b",xlab="number of changepoints",ylab="cost",col=cpcol)
  }else if(output==3){
    plot(x=c(1:x$n),y=sgr,type="l",xlab="t",ylab="Residual")
  }else if(output==4){
    plot(x=c(1:x$n),y=cmrs,type="l",xlab="t",ylab="Cumulative Residual Square")
  }
}