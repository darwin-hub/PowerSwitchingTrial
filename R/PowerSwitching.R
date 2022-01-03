LogRankTestPowerMedian <- function(m1,m2,n,min,max,reps,alpha,r){
  p=vector()
  for(x in 1:reps){
    rate1 <- log(2)/m1
    rate2 <- log(2)/m2
    t1 <- rexp(n, rate=rate1)
    t2 <- rexp(n*r, rate=rate2)
    c1 <- runif(n, min, max)
    c2 <- runif(n*r, min, max)
    group <- c(rep(1,n), rep(2,n*r))
    time <- c(ifelse(t1<c1, t1, c1), ifelse(t2<c2, t2, c2))
    status <- c(ifelse(t1<c1, 1, 0), ifelse(t2<c2, 1, 0))
    diff<-survdiff(Surv(time,status)~group)$chisq
    pvalue<-pchisq(diff,1,lower.tail = FALSE)
    p=c(p,pvalue)
  }
  mean(p<alpha)
}
LogRankTestMedian <- function(m1,m2,min,max,reps,alpha,r,lower,upper,power){
  pwLower <- LogRankTestPowerMedian(m1,m2,n=lower,min,max,reps,alpha,r)
  pwUpper <- LogRankTestPowerMedian(m1,m2,n=upper,min,max,reps,alpha,r)
  K <- 1
  repeat {
    middle <- ceiling((upper+lower)/2)
    pwMiddle <- LogRankTestPowerMedian(m1,m2,n=middle,min,max,reps,alpha,r)
    if ((pwMiddle >= power & middle==lower) | (pwMiddle >= power & middle==upper) | K >= 50) {
      print(c(N=upper, K=K))
      break
    } else if((pwLower > power & pwUpper > power) | (pwLower < power & pwUpper < power)){
      print(rbind(pwLower,pwUpper,lower,upper))
      break
    } else if(pwMiddle >= power & middle!=lower & middle!=upper) {
      upper <- middle; pwUpper <- pwMiddle
    } else if(pwMiddle < power & middle!=upper & pwUpper >= power) {
      lower <- middle; pwLower <- pwMiddle
    }
    K <- K + 1
  }
}
LogRankTestMix2PowerMedian <- function(m1,m2,n,reps,min,max,proportion,s,alpha,r,a,b,random){
  if(random==FALSE){
    p=vector()
    for(x in 1:reps){
      rate1 <- log(2)/m1
      rate2 <- log(2)/m2
      followUp <- sample(c("stay", "leave"), size = n, replace = TRUE, prob = c(1-proportion,proportion))
      tx1 <- rexp(n, rate=rate1)
      c1 <- runif(n, min, max)
      dftx1 <- data.frame(tx1,c1,followUp)
      newDF <- split(dftx1, dftx1$followUp)
      stayT1 <- newDF$stay$tx1
      stayC1 <- newDF$stay$c1
      A <- rate1/rate2
      oldLeaveT1 <- newDF$leave$tx1
      leaveC1 <- newDF$leave$c1
      newLeaveT1 <- c(ifelse(apply(newDF$leave, 1, FUN = min)>s, s+(oldLeaveT1-s)*A, oldLeaveT1))
      timeLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, newLeaveT1, leaveC1))
      statusLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, 1, 0))
      timeStaytx1 <- c(ifelse(stayT1<stayC1, stayT1, stayC1))
      statusStaytx1 <- c(ifelse(stayT1<stayC1, 1, 0))
      newtx2 <- rexp(n*r, rate=rate2)
      newc2 <- runif(n*r, min, max)
      newtimetx2 <- c(ifelse(newtx2<newc2, newtx2, newc2))
      newstatustx2 <- c(ifelse(newtx2<newc2, 1, 0))
      group <- c(rep("mix",n), rep(2,n*r))
      time <- c(c(timeLeavetx1,timeStaytx1),newtimetx2)
      status <- c(c(statusLeavetx1,statusStaytx1),newstatustx2)
      df <- data.frame(group,time,status)
      diff <- survdiff(Surv(time,status)~group, data=df)$chisq
      pvalue<-pchisq(diff,1,lower.tail = FALSE)
      p=c(p,pvalue)
    }
    mean(p<alpha)
  }
  else{
    p=vector()
    for(x in 1:reps){
      rate1 <- log(2)/m1
      rate2 <- log(2)/m2
      followUp <- sample(c("stay", "leave"), size = n, replace = TRUE, prob = c(1-proportion,proportion))
      tx1 <- rexp(n, rate=rate1)
      c1 <- runif(n, min, max)
      dftx1 <- data.frame(tx1,c1,followUp)
      newDF <- split(dftx1, dftx1$followUp)
      stayT1 <- newDF$stay$tx1
      stayC1 <- newDF$stay$c1
      A <- rate1/rate2
      oldLeaveT1 <- newDF$leave$tx1
      leaveC1 <- newDF$leave$c1
      newLeaveT1 <- c(ifelse(apply(newDF$leave, 1, FUN = min)>(time=oldLeaveT1*rbeta(1,a,b)), time+(oldLeaveT1-time)*A, oldLeaveT1))
      timeLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, newLeaveT1, leaveC1))
      statusLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, 1, 0))
      timeStaytx1 <- c(ifelse(stayT1<stayC1, stayT1, stayC1))
      statusStaytx1 <- c(ifelse(stayT1<stayC1, 1, 0))
      newtx2 <- rexp(n*r, rate=rate2)
      newc2 <- runif(n*r, min, max)
      newtimetx2 <- c(ifelse(newtx2<newc2, newtx2, newc2))
      newstatustx2 <- c(ifelse(newtx2<newc2, 1, 0))
      group <- c(rep("mix",n), rep(2,n*r))
      time <- c(c(timeLeavetx1,timeStaytx1),newtimetx2)
      status <- c(c(statusLeavetx1,statusStaytx1),newstatustx2)
      df <- data.frame(group,time,status)
      diff <- survdiff(Surv(time,status)~group, data=df)$chisq
      pvalue<-pchisq(diff,1,lower.tail = FALSE)
      p=c(p,pvalue)
    }
    mean(p<alpha)
  }
}
LogRankTestMix2NMedian <- function(m1,m2,reps,min,max,proportion,s,alpha,r,a,b,random,upper,lower,power){
  pwLower <- LogRankTestMix2PowerMedian(m1,m2,n=lower,reps,min,max,proportion,s,alpha,r,a,b,random)
  pwUpper <- LogRankTestMix2PowerMedian(m1,m2,n=upper,reps,min,max,proportion,s,alpha,r,a,b,random)
  K <- 1
  repeat {
    middle <- ceiling((upper+lower)/2)
    pwMiddle <- LogRankTestMix2PowerMedian(m1,m2,n=middle,reps,min,max,proportion,s,alpha,r,a,b,random)
    if ((pwMiddle >= power & middle==lower) | (pwMiddle >= power & middle==upper) | K >= 50) {
      print(c(N=upper, K=K))
      break
    } else if((pwLower > power & pwUpper > power) | (pwLower < power & pwUpper < power)){
      print(rbind(pwLower, pwUpper,lower,upper))
      break
    } else if(pwMiddle >= power & middle!=lower & middle!=upper) {
      upper <- middle; pwUpper <- pwMiddle
    } else if(pwMiddle < power & middle!=upper & pwUpper >= power) {
      lower <- middle; pwLower <- pwMiddle
    }
    K <- K + 1
  }
}

