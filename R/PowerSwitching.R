LogRankTestPowerMedian <- function(m1,m2,n,Ta,Te,reps,alpha,r){
  p=vector()
  for(x in 1:reps){
    rate1 <- log(2)/m1
    rate2 <- log(2)/m2
    t1 <- rexp(n, rate=rate1)
    t2 <- rexp(n*r, rate=rate2)
    c1 <- c(Te-runif(n, 0, Ta))
    c2 <- c(Te-runif(n*r, 0, Ta))
    group <- c(rep(1,n), rep(2,n*r))
    time <- c(ifelse(t1<c1, t1, c1), ifelse(t2<c2, t2, c2))
    status <- c(ifelse(t1<c1, 1, 0), ifelse(t2<c2, 1, 0))
    diff<-survdiff(Surv(time,status)~group)$chisq
    pvalue<-pchisq(diff,1,lower.tail = FALSE)
    p=c(p,pvalue)
  }
  mean(p<alpha)
}
LogRankTestMedian <- function(m1,m2,Ta,Te,reps,alpha,r,lower,upper,power){
  pwLower <- LogRankTestPowerMedian(m1,m2,n=lower,Ta,Te,reps,alpha,r)
  pwUpper <- LogRankTestPowerMedian(m1,m2,n=upper,Ta,Te,reps,alpha,r)
  K <- 1
  repeat {
    middle <- ceiling((upper+lower)/2)
    pwMiddle <- LogRankTestPowerMedian(m1,m2,n=middle,Ta,Te,reps,alpha,r)
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


LogRankTestMix2PowerMedian <- function(n, m1, m2, r, Ta, Te, proportion, alpha, 
                                       s.dist, pt, rho,
                                       censor.rate=c("AC.only"),
                                       reps) {
  
  if(is.numeric(s.dist)) s <- s.dist else s <- NULL
  
  rate1 <- log(2)/m1
  rate2 <- log(2)/m2
  
  if(is.null(s)) {
    V <- (1/rate1)^2
    mu <- 1/rate1
    
    if(s.dist=="gamma") {
      ab <- sw.gamma(pt, rho, V=V, mu=mu)
      a <- ab[1]; b <- ab[2]
      X.sfun <- function(n) rgamma(n, shape = a, rate = b)
    } else if(s.dist=="beta") {
      ab <- sw.beta(pt, rho, V=V, mu=mu)
      a <- ab[1]; b <- ab[2]
      X.sfun <- function(n) rbeta(n, shape1 = a, shape2 = b)
    } else if(s.dist=="unif") {
      X.sfun <- function(n) runif(n)
    } else if(s.dist=="indepExp") {
      rate.p <- rate1/pt
      X.sfun <- function(n) rexp(n, rate = rate.p)
    }
  }
  
  
  Ta0 <- Ta
  if(Ta==0) Ta <- 0.0001
  
  censoring2.fun <- function(lambda1, Ta, Te) {
    b2 <- (exp(-lambda1*(Te-Ta))-exp(-lambda1*Te))
    b2/(Ta*lambda1)
    #exp(-lambda1*Te)
  }
  h <- NULL
  cens.rate <- censoring2.fun(rate1, Ta, Te)
  
  if(is.numeric(censor.rate)) {
    if(censor.rate <= cens.rate) {
      stop(paste0("censor.rate must be larger than ", round(cens.rate,3)))
    } else {
      censoring.fun <- function(h, lambda1, Ta, Te) {
        if(0 < Te-h & Te-h < Ta) {
          b1 <- (Te-Ta) * exp(-lambda1*(Te-Ta)) - h * exp(-lambda1*h)
          b2 <- (exp(-lambda1*(Te-Ta))-exp(-lambda1*h))
          ( (1 - exp(-lambda1*h)) * (Te-h) + (Ta- (Te-h)) - b2/lambda1)/(Ta*lambda1*h) +
            (h*b2 - (b1 + b2/lambda1))/(Ta*lambda1*h)
        } else if(Ta <= Te-h) {
          1/(h*lambda1) * (1 - exp(-lambda1*h))
        } else {
          b1 <- (Te-Ta) * exp(-lambda1*(Te-Ta)) - Te * exp(-lambda1*Te)
          b2 <- (exp(-lambda1*(Te-Ta))-exp(-lambda1*Te))
          ( Ta - b2/lambda1 + h*b2 - (b1 + b2/lambda1) )/(Ta*lambda1*h)
          #(1 - (1 - h*lambda1 + lambda1*Te) * exp(-lambda1*Te))/(lambda1*h)
        }
      }
      fc <- function(h) censoring.fun(h, rate1, Ta, Te) - censor.rate
      h <- round(uniroot(fc, interval = c(0.001, 10*Te), extendInt = c("yes"))$root, 3)
      cens.rate <- censoring.fun(h, rate1, Ta, Te)
    }
  }
  
  Ta <- Ta0
  
  set.seed(2022)
  p <- vector()
  E1 <- vector()
  E2 <- vector()
  if(!is.null(s) | proportion==0) {

    for(x in 1:reps){
      followUp <- sample(c("stay", "leave"), size = n, replace = TRUE, prob = c(1-proportion,proportion))
      tx1 <- rexp(n, rate=rate1)
      
      if(!is.numeric(censor.rate)) {
        c1 <- c(Te-runif(n, 0, Ta))
      } else {
        v <- runif(n, min=0, max=Ta)
        w1 <- runif(n, min=0, max=h)
        c1 <- ifelse(w1 < Te-v, w1, Te-v)
      }
      
      dftx1 <- data.frame(tx1,c1,followUp)
      newDF <- split(dftx1, dftx1$followUp)
      stayT1 <- newDF$stay$tx1
      stayC1 <- newDF$stay$c1
      A <- rate1/rate2
      oldLeaveT1 <- newDF$leave$tx1
      leaveC1 <- newDF$leave$c1
      if(length(oldLeaveT1)!=0){
        newLeaveT1 <- c(ifelse(apply(newDF$leave[,1:2], 1, FUN = min)>s, s+(oldLeaveT1-s)*A, oldLeaveT1))
      } else{
        newLeaveT1 <- oldLeaveT1
      }
      timeLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, newLeaveT1, leaveC1))
      statusLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, 1, 0))
      timeStaytx1 <- c(ifelse(stayT1<stayC1, stayT1, stayC1))
      statusStaytx1 <- c(ifelse(stayT1<stayC1, 1, 0))
      newtx2 <- rexp(n*r, rate=rate2)
      
      if(!is.numeric(censor.rate)) {
        newc2 <- c(Te-runif(n*r, 0, Ta))
      } else {
        v <- runif(n*r, min=0, max=Ta)
        w1 <- runif(n*r, min=0, max=h)
        newc2 <- ifelse(w1 < Te-v, w1, Te-v)
      }
      
      newtimetx2 <- c(ifelse(newtx2<newc2, newtx2, newc2))
      newstatustx2 <- c(ifelse(newtx2<newc2, 1, 0))
      group <- c(rep("mix",n), rep(2,n*r))
      time <- c(c(timeLeavetx1,timeStaytx1),newtimetx2)
      status <- c(c(statusLeavetx1,statusStaytx1),newstatustx2)
      df <- data.frame(group,time,status)
      diff <- survdiff(Surv(time,status)~group, data=df)$chisq
      pvalue<-pchisq(diff,1,lower.tail = FALSE)
      p=c(p, pvalue)
      E1=c(E1, sum(c(statusLeavetx1,statusStaytx1)))
      E2=c(E2, sum(newstatustx2))
    }
    power <- mean(p<alpha)
    E1 <- mean(E1)
    E2 <- mean(E2)
    
  } else {
    #p=vector()
    for(x in 1:reps){
      followUp <- sample(c("stay", "leave"), size = n, replace = TRUE, prob = c(1-proportion,proportion))
      tx1 <- rexp(n, rate=rate1)

      if(!is.numeric(censor.rate)) {
        c1 <- c(Te-runif(n, 0, Ta))
      } else {
        v <- runif(n, min=0, max=Ta)
        w1 <- runif(n, min=0, max=h)
        c1 <- ifelse(w1 < Te-v, w1, Te-v)
      }
      
      dftx1 <- data.frame(tx1,c1,followUp)
      newDF <- split(dftx1, dftx1$followUp)
      stayT1 <- newDF$stay$tx1
      stayC1 <- newDF$stay$c1
      A <- rate1/rate2
      oldLeaveT1 <- newDF$leave$tx1
      leaveC1 <- newDF$leave$c1
      if(length(oldLeaveT1)!=0){
        if(s.dist=="indepExp") {
          s <- X.sfun(length(oldLeaveT1))
        } else {
          s <- oldLeaveT1 * X.sfun(length(oldLeaveT1))
        }  
        newLeaveT1 <- c(ifelse(apply(newDF$leave[,1:2], 1, FUN = min)>s, s+(oldLeaveT1-s)*A, oldLeaveT1))
      } else {
        newLeaveT1 <- oldLeaveT1
      }
      timeLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, newLeaveT1, leaveC1))
      statusLeavetx1 <- c(ifelse(newLeaveT1<leaveC1, 1, 0))
      timeStaytx1 <- c(ifelse(stayT1<stayC1, stayT1, stayC1))
      statusStaytx1 <- c(ifelse(stayT1<stayC1, 1, 0))
      newtx2 <- rexp(n*r, rate=rate2)
      newc2 <- c(Te-runif(n*r, 0, Ta))
      
      if(!is.numeric(censor.rate)) {
        newc2 <- c(Te-runif(n*r, 0, Ta))
      } else {
        v <- runif(n*r, min=0, max=Ta)
        w1 <- runif(n*r, min=0, max=h)
        newc2 <- ifelse(w1 < Te-v, w1, Te-v)
      }
      
      newtimetx2 <- c(ifelse(newtx2<newc2, newtx2, newc2))
      newstatustx2 <- c(ifelse(newtx2<newc2, 1, 0))
      group <- c(rep("mix",n), rep(2,n*r))
      time <- c(c(timeLeavetx1,timeStaytx1),newtimetx2)
      status <- c(c(statusLeavetx1,statusStaytx1),newstatustx2)
      df <- data.frame(group,time,status)
      diff <- survdiff(Surv(time,status)~group, data=df)$chisq
      pvalue<-pchisq(diff,1,lower.tail = FALSE)
      p=c(p, pvalue)
      E1=c(E1, sum(c(statusLeavetx1,statusStaytx1)))
      E2=c(E2, sum(newstatustx2))
    }
    power <- mean(p<alpha)
    E1 <- mean(E1)
    E2 <- mean(E2)
  }
  
  list(power=power, E1=E1, E2=E2, h=h)
  
}



LogRankTestMix2NMedian <- function(power, lower, upper, m1, m2, r, Ta, Te, proportion, alpha, 
                                   s.dist, pt, rho,
                                   censor.rate=c("AC.only"),
                                   reps) {

  pwLower <- LogRankTestMix2PowerMedian(n=lower, m1, m2, r, Ta, Te, proportion, alpha, 
                                        s.dist, pt, rho,
                                        censor.rate,
                                        reps=100)$power
  pwUpper <- LogRankTestMix2PowerMedian(n=upper, m1, m2, r, Ta, Te, proportion, alpha, 
                                        s.dist, pt, rho,
                                        censor.rate,
                                        reps=100)$power
  K <- 1
  repeat {
    middle <- ceiling((upper+lower)/2)
    Middle <- LogRankTestMix2PowerMedian(n=middle, m1, m2, r, Ta, Te, proportion, alpha, 
                                         s.dist, pt, rho,
                                         censor.rate,
                                         reps=100)
    pwMiddle <- Middle$power
    
    if ((pwMiddle >= power & middle==lower) | (pwMiddle >= power & middle==upper) | K >= 50) {
      N <- upper
      break
    } else if((pwLower > power & pwUpper > power) | (pwLower < power & pwUpper < power)){
      print(rbind(pwLower, pwUpper,lower,upper))
      print("Please increase upper (N).")
      N <- NA
      break
    } else if(pwMiddle >= power & middle!=lower & middle!=upper) {
      upper <- middle; pwUpper <- pwMiddle
    } else if(pwMiddle < power & middle!=upper & pwUpper >= power) {
      lower <- middle; pwLower <- pwMiddle
    }
    K <- K + 1
  }

  
  if (!is.na(N)) {
    #lower <- N-20; upper <- N+20
    lower <- 0.9*N; upper <- 1.1*N
    
    K <- 1
    repeat {
      middle <- ceiling((upper+lower)/2)
      Middle <- LogRankTestMix2PowerMedian(n=middle, m1, m2, r, Ta, Te, proportion, alpha, 
                                           s.dist, pt, rho,
                                           censor.rate,
                                           reps)
      pwMiddle <- Middle$power
      
      if ((pwMiddle >= power & middle==lower) | (pwMiddle >= power & middle==upper) | K >= 50) {
        N=upper; K=K; E1=Middle$E1; E2=Middle$E2; h=Middle$h
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
    c(N1=N, N2=r*N, E1=E1, E2=E2, K=K, h=h)
  }
  
  
}




sw.gamma <- function(p, rho, V, mu) {

  if (p<0 | abs(rho)>1) {
    stop(print("pt must be larger than 0, and rho must be between -1 and 1"))
  } else {
    fb <- function(b) {
      a <- p * b
      mu.x <- p
      V.x <- a/b^2
      num <- mu.x * sqrt(V)
      den <- sqrt((V.x + mu.x^2) * V + V.x * mu^2)
      num/den
    }
    min.rho <- min(fb(c(seq(0.001, 0.01, 0.001))))
    
    if (rho < min.rho) {
      print(paste0("rho must be larger than ", round(min.rho,2), " under p = ", p))
    } else {
      fb.rho <- function(b) fb(b) - rho
      
      if(abs(fb.rho(0.01))<0.001) {
        b <- 0.01
        a <- p*b
      } else {
        solve.b <- uniroot(fb.rho, interval = c(0.01, 20), extendInt="yes")
        b <- solve.b$root;
        a <- p*b
      }  
      c(a=a, b=b)
    }  
  }
  
}



sw.beta <- function(p, rho, V, mu) {

  if (p<0 | p>1 | abs(rho)>1) {
    stop(print("pt must be between 0 and 1, and rho must be between -1 and 1"))
  } else {
    fb <- function(b) {
      a <- p/(1-p)*b
      mu.x <- p
      V.x <- a*b/((a+b)^2 * (a+b+1))
      num <- mu.x * sqrt(V)
      den <- sqrt((V.x + mu.x^2) * V + V.x * mu^2)
      num/den #- rho
    }
    min.rho <- min(fb(c(seq(0.001, 0.01, 0.001))))
    
    if (rho < min.rho) {
      stop(print(paste0("rho must be larger than ", round(min.rho,2), " under pt = ", p)))
    } else {
      fb.rho <- function(b) fb(b) - rho
      
      if(abs(fb.rho(0.01))<0.001) {
        b <- 0.01
        a <- p/(1-p)*b
      } else {
        solve.b <- uniroot(fb.rho, interval = c(0.01, 20), extendInt="yes")
        b <- solve.b$root;
        a <- p/(1-p)*b
      }  
      c(a=a, b=b)
    }  
  }
  
}
