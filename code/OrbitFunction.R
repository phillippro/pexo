#update encounter model
#########################
#inf<- 1000
library(cubature)
library(HyperbolicDist)
library(base)
library(pracma)
library(ks)
library(numDeriv)
mypeaks1<- function(x,span){
       index.min<- as.logical(mat.or.vec(length(x),1))
       index.max<- as.logical(mat.or.vec(length(x),1))
       j<- round(span/2)
       d<- 10*(max(x)-min(x))/length(x)
        for(k in 1:length(x)){
          if(any(k==1:j) || any(k==(length(x)-j+1):length(x))){
          index.min[k]<- FALSE
          index.max[k]<- FALSE
          } else
            if(all(x[k]>x[(k-j):(k-1)],x[k]>x[(k+1):(k+j)],abs(x[(k-j+1):k]-x[(k-j):(k-1)])<d,abs(x[k:(k+j-1)]-x[(k+1):(k+j)])<d)){
              index.min[k]<-FALSE
              index.max[k]<-TRUE
            } else
              if(all(x[k]<x[(k-j):(k-1)],x[k]<x[(k+1):(k+j)],abs(x[(k-j+1):k]-x[(k-j):(k-1)])<d,abs(x[k:(k+j-1)]-x[(k+1):(k+j)])<d)){
              index.min[k]<-TRUE
              index.max[k]<-FALSE
              } else{
                index.min[k]<- FALSE
                index.max[k]<- FALSE
}}
     return(list(index.min=index.min,index.max=index.max))
}

mypeaks2<- function(x,span){
       index.min<- as.logical(mat.or.vec(length(x),1))
       index.max<- as.logical(mat.or.vec(length(x),1))
       j<- round(span/2)
       d<- 10*(max(x)-min(x))/length(x)
       for(k in 1:length(x)){
          if(any(k==1:j) || any(k==(length(x)-j+1):length(x))){
          index.min[k]<- FALSE
          index.max[k]<- FALSE
          } else
            if(all(x[(k-j+1):k]>=x[(k-j):(k-1)],x[k:(k+j-1)]>=x[(k+1):(k+j)])){
              index.min[k]<-FALSE
              index.max[k]<-TRUE
            } else
              if(all(x[(k-j+1):k]<=x[(k-j):(k-1)],x[k:(k+j-1)]<=x[(k+1):(k+j)])){
              index.min[k]<-TRUE
              index.max[k]<-FALSE
              } else{
                index.min[k]<- FALSE
                index.max[k]<- FALSE
}}
     return(list(index.min=index.min,index.max=index.max))
}

mypeaks3<- function(x,span){
       index.min<- as.logical(mat.or.vec(length(x),1))
       index.max<- as.logical(mat.or.vec(length(x),1))
       j<- round(span/2)
       d<- 10*(max(x)-min(x))/length(x)
       for(k in 1:length(x)){
          if(any(k==1:j) || any(k==(length(x)-j+1):length(x))){
          index.min[k]<- FALSE
          index.max[k]<- FALSE
          } else
            if(all(x[(k-j+1):k]>x[(k-j):(k-1)],x[k:(k+j-1)]>x[(k+1):(k+j)])){
              index.min[k]<-FALSE
              index.max[k]<-TRUE
            } else
              if(all(x[(k-j+1):k]<x[(k-j):(k-1)],x[k:(k+j-1)]<x[(k+1):(k+j)])){
              index.min[k]<-TRUE
              index.max[k]<-FALSE
              } else{
                index.min[k]<- FALSE
                index.max[k]<- FALSE
}}
     return(list(index.min=index.min,index.max=index.max))
}

################################################
##########Define the morphology of massless arm
################################################
radiu.start<- function(theta,par){
            Rmin<- par$Rmin
            thetamin<- par$thetamin
            alpha<- par$alpha
            dR<- par$dR
            rvalue<- Rmin*exp((theta-thetamin)/alpha)+dR/2
            return(rvalue)
}
radiu.center<- function(theta,par){
            Rmin<- par$Rmin
            thetamin<- par$thetamin
            alpha<- par$alpha
            rvalue<-Rmin*exp((theta-thetamin)/alpha)
            return(rvalue)
}
radiu.end<- function(theta,par){
            Rmin<- par$Rmin
            thetamin<- par$thetamin
            alpha<- par$alpha
            dR<- par$dR
            rvalue<- Rmin*exp((theta-thetamin)/alpha)-dR/2
            return(rvalue)
}

####reverse function
phi.start <- function(R,par){
            Rmin<- par$Rmin
            thetamin<- par$thetamin
            alpha<- par$alpha
            dR<- par$dR
           theta<- alpha*log((R-dR/2)/Rmin)+thetamin
           return(theta)
}
phi.center <- function(R,par){
            Rmin<- par$Rmin
            thetamin<- par$thetamin
            alpha<- par$alpha
           theta<- alpha*log(R/Rmin)+thetamin
           return(theta)
}
phi.end <- function(R,par){
            Rmin<- par$Rmin
            thetamin<- par$thetamin
            alpha<- par$alpha
            dR<- par$dR
           theta<- alpha*log((R+dR/2)/Rmin)+thetamin
           return(theta)
}

#draw a circle: the bulge
circle<- function(R){
         phi<-seq(0,2*pi,0.01)
         x<- R*cos(phi)
         y<- R*sin(phi)
         return(list(x=x,y=y))
}
###########################################
###count the crossings of a arm
###########################################
#judge a crossing events
index.arm <- function(R0,phi0,R1,phi1,par){
              Rmin<- par$Rmin
              Rmax<- par$Rmin+par$extent
              if(R0>Rmax || R0<Rmin || R1>Rmax || R1<Rmin){
                   ind.start<- FALSE
                   ind.center<- FALSE
                   ind.end<- FALSE} else
                   if(R0>radiu.start(phi0,par) & R1<=radiu.start(phi1,par)){
                       ind.start<- TRUE
                       ind.center<- FALSE
                       ind.end<- FALSE}else
                       if(R0>radiu.center(phi0,par) & R1<=radiu.center(phi1,par)){
                           ind.start<- FALSE
                           ind.center<- TRUE
                           ind.end<- FALSE} else
                           if(R0>radiu.end(phi0,par) & R1<=radiu.end(phi1,par)){
                               ind.start<- FALSE
                               ind.center<- FALSE
                               ind.end<- TRUE} else{
                                                    ind.start<- FALSE
                                                    ind.center<- FALSE
                                                    ind.end<- FALSE}
               return(list(start=ind.start,center=ind.center,end=ind.end))
}

#find the index of the crossing points for a specific arm model
index.cross<- function(R,phi,par){
    thetamin <- par$thetamin
    Rmax <- par$Rmin+par$extent
    thetamax <- phi.center(Rmax,par)
    index.start <- as.logical(mat.or.vec(length(phi),1))
    index.center <- as.logical(mat.or.vec(length(phi),1))
    index.end <- as.logical(mat.or.vec(length(phi),1))
    phi <- phi%%(2*pi)
########judge the crossing events
    for(k in 1:(length(phi)-1)){
        if(phi[k]>thetamin & phi[k]<thetamax & (phi[k]+2*pi)>thetamax){
            index.start[k] <- index.arm(R[k],phi[k],R[k+1],phi[k+1],par)$start
            index.center[k] <- index.arm(R[k],phi[k],R[k+1],phi[k+1],par)$center
            index.end[k] <- index.arm(R[k],phi[k],R[k+1],phi[k+1],par)$end}else
                                                                              if(phi[k]+2*pi<=thetamax){
                                                                                  index.start[k] <- index.arm(R[k],phi[k]+2*pi,R[k+1],phi[k+1]+2*pi,par)$start
                                                                                  index.center[k] <- index.arm(R[k],phi[k]+2*pi,R[k+1],phi[k+1]+2*pi,par)$center
                                                                                  index.end[k] <- index.arm(R[k],phi[k]+2*pi,R[k+1],phi[k+1]+2*pi,par)$end
                                                                              }else{
                                                                                  index.start[k]<- FALSE
                                                                                  index.center[k]<- FALSE
                                                                                  index.end[k]<- FALSE
                                                                              }
    }
############The last position of the Sun should not be the crossing event
    index.start[length(phi)]<- FALSE
    index.center[length(phi)]<- FALSE
    index.end[length(phi)]<- FALSE
    return(list(start=index.start,center=index.center,end=index.end))
}

#In order to fill the lossing crossing events, we have to make a sequence of
#crossing time according to t1->t2->t3->t1->t2->t3->t1->t2->t3->.... If some
#crossing part (start, center, end) are lost, e.g. t2(center crossing time) is lost, then we could code to
#fill it to be NA. Then, we know which crossing part is lost. Based on this
#work, we could know how which arm is lost. So there are two steps to do:
#1. find the losing crossing part (start, center, end)
#2. find the losing crossing arm (arm1, arm2, (arm3, arm4... for 4-arm model))
#Find the lossing crossing part
order.cross <- function(time.start,time.center,time.end){
            time <- sort(c(time.start,time.center,time.end))
            i <- 1
            while(i<=(length(time)-2)){
              if(is.element(time[i],time.start)){
                  if(is.element(time[i+1],time.center)){
                       if(is.element(time[i+2],time.end)){
                           time <- time
                       }else{
                           time <- append(time, NA, after=i+1)}
                  }else if(is.element(time[i+1], time.end)){
                       time <- append(time, NA, after=i)}else{
                       time <- append(time, rep(NA,2), after=i)}
              }else if(is.element(time[i],time.center)){
                  if(is.element(time[i+1],time.end)){
                  time <- append(time, NA, after=i-1)}else{
                  time <- append(time, NA, after=i-1)
                  time <- append(time, NA, after=i+1)}
              }else{
              time <- append(time, rep(NA,2), after=i-1)
              }
            i <- i+3
            }
            if((length(time)-i)==-1){
              time <- time}else
              if(length(time)-i==0){
                  if(is.element(time[i],time.start)){
                     time <- append(time, rep(NA,2), after=i)}else
                     if(is.element(time[i],time.center)){
                         time <- append(time, NA, after=i-1)
                         time <- append(time, NA, after=i+1)}else{
                              time <- append(time, rep(NA,2), after=i-1)}
                }else{
                  if(!any(is.element(time[i:length(time)],time.start))){
                      time <- append(time, NA, after=i-1)}else
                      if(!any(is.element(time[i:length(time)],time.center))){
                         time <- append(time, NA, after=i)}else{
                               time <- append(time, NA, after=i+1)}
                }
           # print(time)
            if(length(time)==0){
                 time.start <- NA
                 time.center <- NA
                 time.end <- NA}else{
            time.start <- time[seq(1,length(time),3)]
            time.center <- time[seq(2,length(time),3)]
            time.end <- time[seq(3,length(time),3)]}
            return(list(start=time.start,center=time.center,end=time.end))
}
#################denoted the missing crossing parts to be NA for variables
var.na<- function(var,standard){
    j<- 1
    while(j <= length(standard)){
       if(is.na(standard[j]) & !is.na(var[j])){
          var<- append(var,NA,after=j-1)}else{
          var<- var}
       j<- j+1}
 return(var)
}
#################denoted the missing arm crossings to be NA, the default
#sequence is 2->1->2->1->2->...
label.arm <- function(label.start, label.center, label.end){
             lab <- rep(NA, length(label.start))
             for(i in 1:length(lab)){
                 temp <- c(label.start[i],label.center[i],label.end[i])
                 if(all(is.na(temp))){
                    lab[i] <- NA}else{
                           lab[i] <- temp[!is.na(temp)][1]}
             }

             if(!is.na(lab[1]) & lab[1]==1){
             lab <- c(NA, lab)}
             #print(lab)
             j <- 1
             while(j <= (length(lab)-1)){
                 if(!is.na(lab[j]) & !is.na(lab[j+1]) & lab[j]==lab[j+1]){
                          label.start <- append(label.start, NA, after=j)
                          label.center <- append(label.center, NA, after=j)
                          label.end <- append(label.end, NA, after=j)
                          lab <- append(lab, NA, after=j)
                         # print(list(label.start,label.center,label.end,lab))
                          }
             j <- j+1
             }

             return(list(start=label.start, center=label.center, end=label.end))
}

# Given unnormalized TS model function func with parameters modpar, calculate and return
# normalization constant for function extending over time span given by variable time.
# The normalization constant is what we multiply the function by such that it's integral is unity,
# so we return the *inverse* of the integral of the function.
# If the integral is zero, then return zero (rather than infinity). This indicates that the function is non-normalizable.
# Note that I test for 1/val==Inf, not val==0, as the former can be statisfied when the latter is not!
# time=list(min, max)
# Nsamp=number of samples for time calculation
#func:discreate function
CalcNormVal <- function(func, time) {
  val <- as.numeric(1/2*(func[-length(func)]+func[-1])%*%diff(time))
  if(val<0) stop("Error in CalcNormVal: integral is negative\n")
  if(1/val==Inf) { return(0) } else { return(func/val) }
}

#########
# Convolve the normalized function func(t) with a 1D Gaussian N_t(x,xsig) over the time range given by time.
# func = func(t, funcpar), where t is continuous time variable and funcpar are other parameters
# funcpar$normval should be set in advance (is used within func to return normalized probability)
# x = mean of Gaussian (scalar)
# xsig = sd of Gaussian (scalar) (must be positive)
# Approximate the integral with a sum with uniform sampling.
# As the Gaussian will generally be much narrower than the time span, I just sample over that.
# 2011-01-14: I double checked that the integration, esp. the time ranges used and divided by, are and were correct:
# If we now integrated instead over the full time range by doubling the time span and doubling the number of samples (to keep the same
# sampling as in the first case), then we would now add lots of zeros in the sum, leaving it unchanged, and the sample
# width is unchanged because both the time span and the number of samples have doubled. If we just doubled the time span but
# kept Nsamp the same, then the sum would be half, because we now have half as many samples, so again integ is the same.
ConvGauss.norm <- function(x, xsig, func, normval,time) {
  Nsig <- 6              # assume Gaussian is zero for |t-x|/xsig > Nsig
  Nsamp <- 2*10*Nsig + 1 # number of samples for integral (here set to well sample range ±Nsig)
  t.lower <- max(min(time), x-Nsig*xsig)
  t.upper <- min(max(time), x+Nsig*xsig)
  if(t.lower>=max(time) || t.upper<=min(time)) {return(0)}
  t <- seq(from=t.lower, to=t.upper, length.out=Nsamp)
  integ <- sum(normval*func(t)*dnorm(t, mean=x, sd=xsig,log=FALSE))# sum over product of two normalized functions
  integ <- integ*(t.upper-t.lower)/Nsamp  # (t.upper-t.lower)/Nsamp is the width of each time sample in sum
  return(integ)
}

ConvGauss <- function(x, xsig, func, time) {
  Nsig <- 6              # assume Gaussian is zero for |t-x|/xsig > Nsig
  Nsamp <- 2*10*Nsig + 1 # number of samples for integral (here set to well sample range ±Nsig)
  t.lower <- max(min(time), x-Nsig*xsig)
  t.upper <- min(max(time), x+Nsig*xsig)
  if(t.lower>=max(time) || t.upper<=min(time)) {return(0)}
  t <- seq(from=t.lower, to=t.upper, length.out=Nsamp)
  integ <- sum(func(t)*dnorm(t, mean=x, sd=xsig,log=FALSE))# sum over product of two normalized functions
  integ <- integ*(t.upper-t.lower)/Nsamp  # (t.upper-t.lower)/Nsamp is the width of each time sample in sum
  return(integ)
}

ConvFunc <- function(func1,func2,t){
  integ <- sum(func1(t)*func2(t))# sum over product of two normalized functions
  integ <- integ*(max(t)-min(t))/length(t)  # (t.upper-t.lower)/Nsamp is the width of each time sample in sum
  return(integ)
}

# Convolve the normalized function func(t) with a 1D Uniform distribution extending from x1:x2 over the time range given by time.
# funcpar$normval should be set in advance (is used within func to return normalized probability)
# Approximate the integral with a sum with uniform sampling at Nsamp points
ConvUniform <- function(t.cen, dur, func, time) {
  if(t.cen<=0) stop("Error in ConvUniform: x must be positive\n")
  Nsamp <- as.integer(5*dur)  # appropriate if time is in Myr
  x1<- t.cen-dur/2
  x2<- t.cen+dur/2
  t.lower <- max(min(time), x1)
  t.upper <- min(max(time), x2)
  t <- seq(from=x1, to=x2, length.out=Nsamp)
  ## In full:
  #  integ <- (1/dur) * sum( func(t, funcpar) )  # sum over product of two normalized functions
  #  integ <- integ*x/Nsamp                    # x/Nsamp is width of each time sample in sum
  integ <- sum( func(t) ) / Nsamp
  return(integ)
}

#Uniform Model
# Constant probability
# parameters: normval
uni.prob <- function(t, modpar) {
  return(rep(modpar$normval, length(t)))
}

# Given unnormalized TS model function func with parameters modpar, calculate and return
# normalization constant for function extending over time span given by variable time.
# The normalization constant is what we multiply the function by such that it's integral is unity,
# so we return the *inverse* of the integral of the function.
# If the integral is zero, then return zero (rather than infinity). This indicates that the function is non-normalizable.
# Note that I test for 1/val==Inf, not val==0, as the former can be statisfied when the latter is not!
# time=list(min, max)
# Nsamp=number of samples for time calculation
####calculate the contineous function
CalcNormVal.cont <- function(func, time, Nsamp) {
  t <- seq(from=time$min, to=time$max, length.out=Nsamp)
  val <- sum(func(t))
  val <- (time$max-time$min)*val/Nsamp
  if(val<0) stop("Error in CalcNormVal: integral is negative\n")
  if(1/val==Inf) { return(0) } else { return(1/val) }
}

simEnc <- function(tmax,dmax,rs,vs,rate=20,vary.apex=FALSE){
###dmax [kpc]; tmax [Myr]; rs [kpc]; vs [kpc/Myr]
    kpcmyr2auyr <- 1e3*206265/1e6
    kpc2au <- 1e3*206265
    auyr2kms <- 4.74047
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    ###peculiar velocity
    Rs <- sqrt(rs[1]^2+rs[2]^2)
    vc0 <- -Vcirc(Rs,par.gal)#kpc/myr
    phi <- xy2phi(rs[1],rs[2])
    Vc <- vc0*c(-sin(phi),cos(phi),0)
    vp0 <- (vs-Vc)*kpcmyr2kms
###
    vp.enc <- sqrt(vp0%*%vp0)
    rate <- round(rate*(dmax/0.001)^2*(vp.enc/vp.sun))
    if(FALSE){
        cat('vp.enc=',vp.enc,'km/s\n')
        cat('vp.sun=',vp.sun,'km/s\n')
        cat('Rs=',Rs,'kpc\n')
        cat('Vc=',Vc*kpcmyr2kms,'km/s\n')
        cat('vs=',vs*kpcmyr2kms,'km/s\n')
        cat('vp=',vp0,'km/s\n')
        cat('rate=',rate,'\n\n')
    }
    if(tmax<1){
        tenc <- runif(rate,0,1)
        ind <- which(tenc<tmax)
        tenc <- tenc[ind]
    }else{
        tenc <- runif(rate*tmax,0,tmax)
    }
    Nenc <- length(tenc)
    de <- runif(Nenc,0,1)^0.5*dmax
    Ntry <- 10*Nenc
    type.name <- c("B0","A0","A5","F0","F5","G0","G5","K0","K5","M0","M5","wd","gi")
    mass <- c(9,3.2,2.1,1.7,1.3,1.1,0.93,0.78,0.69,0.47,0.21,0.9,4)
    enc.freq <- c(0.005,0.03,0.04,0.15,0.08,0.22,0.35,0.34,0.85,1.29,6.39,0.72,0.06)
    vsig.apex <- c(14.7,19.7,23.7,29.1,36.2,37.4,39.2,34.1,43.4,42.7,41.8,63.4,41.0)
    vs.apex <- c(18.6,17.1,13.7,17.1,17.1,26.4,23.9,19.8,25.0,17.3,23.3,38.3,21.0)
    flags <- c()
    for(j in 1:length(type.name)){
        flags <- c(flags,rep(j,round(enc.freq[j]*1e3)))
    }
#####initial tries
    ii <- sample(flags,Ntry,replace=TRUE)
    Menc <- mass[ii]
    vpp <- sqrt(sum(vp0^2))
    if(vary.apex){
        vp <- t(outer(vp0,vs.apex[ii]))/sqrt(sum(vp.sun^2))
    }else{
        vp <- vp0
    }
    eta <- rnorm(length(ii),0,1)
    nu <- rnorm(length(ii),0,1)
    w <- rnorm(length(ii),0,1)
    vstar <- sqrt(1/3)*cbind(vsig.apex[ii]*eta,vsig.apex[ii]*nu,vsig.apex[ii]*w)
    if(!is.null(dim(vp))){
        venc <- vstar-vp
    }else{
        venc <- t(t(vstar)-vp)
    }
    ve <- sqrt(rowSums(venc^2))#mode
    V <- vpp+3*vsig.apex[ii]
    beta <- runif(Ntry,0,1)
    ind <- which((beta*V)<ve)
    index <- ind[1:Nenc]
####encounters following P(venc)~venc
    ve <- ve[index]
    venc <- venc[index,]
    Menc <- Menc[index]
    blv <- xyz2bl.vec(venc[,1],venc[,2],venc[,3])
    b.enc <- blv[,1]
    l.enc <- blv[,2]
    a.enc <- runif(Nenc,0,2*pi)
    renc <- s2g(0,de*cos(a.enc),de*sin(a.enc),-b.enc,-l.enc)
    re <- sqrt(rowSums(renc^2))
    if(vary.apex){
        vv <- vp[index,]
        cosk <- rowSums(renc*vv)/(re*sqrt(rowSums(vv*vv)))
        cosa <- rowSums(venc*vv)/(ve*sqrt(rowSums(vv*vv)))
    }else{
        cosk <- renc%*%vp/(re*vpp)
        cosa <- venc%*%vp/(ve*vpp)
    }
    return(data.frame(tenc=tenc,x=renc[,1],y=renc[,2],z=renc[,3],vx=venc[,1]/kpcmyr2kms,vy=venc[,2]/kpcmyr2kms,vz=venc[,3]/kpcmyr2kms,Menc=Menc,cosk=cosk,cosa=cosa,ve=ve,re=re))
}

simEnc.fixFi <- function(tmax,dmax,rs,vs,rate=20,vary.apex=FALSE){
###dmax [kpc]; tmax [Myr]; rs [kpc]; vs [kpc/Myr]
    rate0 <- rate
    kpcmyr2auyr <- 1e3*206265/1e6
    kpc2au <- 1e3*206265
    auyr2kms <- 4.74047
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
###peculiar velocity
    Rs <- sqrt(rs[1]^2+rs[2]^2)
    vc0 <- -Vcirc(Rs,par.gal)#kpc/myr
    phi <- xy2phi(rs[1],rs[2])
    Vc <- vc0*c(-sin(phi),cos(phi),0)
    vp0 <- (vs-Vc)*kpcmyr2kms
###
    vp.enc <- sqrt(vp0%*%vp0)
###different types
    type.name <- c("B0","A0","A5","F0","F5","G0","G5","K0","K5","M0","M5","wd","gi")
    mass <- c(9,3.2,2.1,1.7,1.3,1.1,0.93,0.78,0.69,0.47,0.21,0.9,4)
    enc.freq <- c(0.005,0.03,0.04,0.15,0.08,0.22,0.35,0.34,0.85,1.29,6.39,0.72,0.06)
    vsig.apex <- c(14.7,19.7,23.7,29.1,36.2,37.4,39.2,34.1,43.4,42.7,41.8,63.4,41.0)
    vs.apex <- c(18.6,17.1,13.7,17.1,17.1,26.4,23.9,19.8,25.0,17.3,23.3,38.3,21.0)
###
#    rate <- round(rate*(dmax/0.001)^2*vp.enc/vp.sun)
    rate <- round(rate0*(dmax/0.001)^2*10)#define a hugh encounter rate and select them later
    if(FALSE){
        cat('dmax=',dmax,'\n')
        cat('vp.enc=',vp.enc,'km/s\n')
        cat('vp.sun=',vp.sun,'km/s\n')
        cat('Rs=',Rs,'kpc\n')
        cat('Vc=',Vc*kpcmyr2kms,'km/s\n')
        cat('vs=',vs*kpcmyr2kms,'km/s\n')
        cat('rate=',rate,'\n\n')
        cat('tmax=',tmax,'\n')
        cat('rate=',rate,'\n')
    }
    if(tmax<1){
        tenc <- runif(rate,0,1)
        ind <- which(tenc<tmax)
        tenc <- tenc[ind]
    }else{
        tenc <- runif(rate*tmax,0,tmax)
    }
    Nenc <- length(tenc)
    de <- runif(Nenc,0,1)^0.5*dmax
    Ntry <- 100
    flags <- c()
    for(j in 1:length(type.name)){
        flags <- c(flags,rep(j,round(enc.freq[j]*1e3)))
    }
#####fix encounter rate for each category
    cosa <- Renc <- Venc <- Menc <- cosk <- Ve <- Re <- c()
    for(k in 1:Nenc){
        ii <- sample(flags,1)
        if(vary.apex){
            vp <- vs.apex[ii]*vp0/sqrt(vp.sun%*%vp.sun)
        }else{
            vp <- vp0
        }
        for(j in 1:Ntry){
            eta <- rnorm(3,0,1)
            vstar <- sqrt(1/3)*eta*vsig.apex[ii]
            if(is.null(dim(vp))){
                vpp <- vp
            }else{
                vpp <- vp[ii,]
            }
            venc <- vstar-vpp
            ve <- sqrt(venc%*%venc)
            V <- sqrt(vpp%*%vpp)+3*vsig.apex[ii]
            vtry <- runif(1,0,V)
            if(vtry<ve) break()
        }
        if(j==Ntry) cat('Maximum try times reached!\n')
        Venc <- rbind(Venc,venc)
        Ve <- rbind(Ve,ve)
        Menc <- c(Menc,mass[ii])
        blv <- xyz2bl.vec(venc[1],venc[2],venc[3])
        b.enc <- blv[1]
        l.enc <- blv[2]
        a.enc <- runif(1,0,2*pi)
        renc <- as.numeric(s2g(0,de[k]*cos(a.enc),de[k]*sin(a.enc),-b.enc,-l.enc))
        re <- sqrt(renc%*%renc)
        Re <- c(Re,re)
        Renc <- rbind(Renc,renc)
        cosk <- c(cosk,renc%*%vp/(sqrt(renc%*%renc)*sqrt(vp%*%vp)))
        cosa <- c(cosa,venc%*%vp/(sqrt(venc%*%venc)*sqrt(vp%*%vp)))
    }
###redefine encounter rate based on simulated encounters; assume 50 km/s is the mean encounter velocity for the Sun
    rate <- round(rate0*(dmax/0.001)^2*mean(Ve)/50)
    index <- sample(1:nrow(Renc),rate)
    return(data.frame(tenc=tenc[index],x=Renc[index,1],y=Renc[index,2],z=Renc[index,3],vx=Venc[index,1]/kpcmyr2kms,vy=Venc[index,2]/kpcmyr2kms,vz=Venc[index,3]/kpcmyr2kms,Menc=Menc[index],cosk=cosk[index],cosa=cosa[index],ve=Ve[index],re=Re[index]))
}
######Uniform model, we can convolve
# Return event epochs within specified time window for a fixed number of events with epochs drawn
# independently from a uniform random distribution. The events are not sorted into increasing time order
# fpar=list(Nevent)
RanUnifIndEvents <- function(vpar=NULL,fpar,time) {
  epochs  <- runif(fpar$Nevent, min=time$min, max=time$max)
  return(epochs)
}

######Gaussian function
##generate several gaussian functions with different means and uncertainties
##integrate them into one function
gaus <- function(x,mean,sigma){
                  1/(sqrt(2*pi)*sigma)*exp(-(x-mean)^2/2/sigma^2)}

gaus.unnorm<- function(x,a,mu,sig){
                  a*exp(-(x-mu)^2/(2*sig^2))}


randb.mean <- function(x,mean,sigma,fracb){
          val<- fracb/((1-fracb)*sqrt(2*pi)*mean(sigma))
          for(i in 1:length(mean)){
                 val<- val+dnorm(x,mean[i],sigma[i])}
          return(val)
}

randb.peak <- function(x,mean,sigma,fracb){
          val<- fracb/((1-fracb)*sqrt(2*pi)*min(sigma))
          for(i in 1:length(mean)){
                 val<- val+dnorm(x,mean[i],sigma[i])}
          return(val)
}

norm.like <- function(x,mean,sigma,b,a){
          val <- b
          for(i in 1:length(mean)){
                 val<- val+a[m]*sqrt(2*pi)*sigma[i]*dnorm(x,mean[i],sigma[i])}
          return(val)
}

Orbit <- function(t,state,parameters){
    with(as.list(c(state,parameters)),{
                                        #rate of change
        dR<- Rdot
        dRdot<- -numerical.derivative(function(R,z) phit(R,z,parameters),R,z)+R*phidot^2
        dphi<- phidot
        dphidot<- -2*Rdot/R*phidot
        dz<- zdot
        dzdot<- -numerical.derivative(function(z,R) phit(R,z,parameters),z,R)
                                        #return the rate of change
        list(c(dR,dRdot,dphi,dphidot,dz,dzdot))
    })
}

Orbit.point <- function(t,state,parameters){
        out <- rep(NA,6*Nt)
        ind <- 1:Nt
        xs <- state[(1:Nt)*6-5]
        ys <- state[(1:Nt)*6-4]
        zs <- state[(1:Nt)*6-3]
        kpc2au <- 1e3*206265
        for(j in 1:length(Ms)){
            out[(j-1)*6+1:3] <- state[(j-1)*6+4:6]
            d3 <- ((xs[-j]-xs[j])^2+(ys[-j]-ys[j])^2+(zs[-j]-zs[j])^2)^1.5
            out[(j-1)*6+4] <- -mu*sum(Ms[-j]*(xs[j]-xs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6
            out[(j-1)*6+5] <- -mu*sum(Ms[-j]*(ys[j]-ys[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6
            out[(j-1)*6+6] <- -mu*sum(Ms[-j]*(zs[j]-zs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6
        }
        list(out)
}

Orbit2.point <- function(t,y){
        out <- rep(NA,6*Nt)
        mu <- (2*pi)^2
        xs <- y[(1:Nt)*6-5]
        ys <- y[(1:Nt)*6-4]
        zs <- y[(1:Nt)*6-3]
        kpc2au <- 1e3*206265
        for(j in 1:length(Ms)){
            out[(j-1)*6+1:3] <- y[(j-1)*6+4:6]
            d3 <- ((xs[-j]-xs[j])^2+(ys[-j]-ys[j])^2+(zs[-j]-zs[j])^2)^1.5
            out[(j-1)*6+4] <- -mu*sum(Ms[-j]*(xs[j]-xs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6
            out[(j-1)*6+5] <- -mu*sum(Ms[-j]*(ys[j]-ys[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6
            out[(j-1)*6+6] <- -mu*sum(Ms[-j]*(zs[j]-zs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6
        }
        out
}

Orbit.tide <- function(t,state,parameters){
    Nt <- length(state)/6
    out <- rep(NA,6*Nt)
    ind <- 1:Nt
    xs <- state[(1:Nt)*6-5]
    ys <- state[(1:Nt)*6-4]
    zs <- state[(1:Nt)*6-3]
    for(j in 1:Nt){
        out[(j-1)*6+1:3] <- state[(j-1)*6+4:6]
        R <- sqrt(sum(state[(j-1)*6+1:2]^2))
        z <- state[(j-1)*6+3]
        phi <- xy2phi(state[(j-1)*6+1],state[(j-1)*6+2])
        ar <- -numerical.derivative(function(R,z) phit(R,z,parameters),R,z)
        aphi <- 0
        az <- -numerical.derivative(function(z,R) phit(R,z,parameters),z,R)
        d3 <- ((xs[-j]-xs[j])^2+(ys[-j]-ys[j])^2+(zs[-j]-zs[j])^2)^1.5
        out[(j-1)*6+4] <- ar*cos(phi)#kpc/Myr^2
        out[(j-1)*6+5] <- ar*sin(phi)
        out[(j-1)*6+6] <- az
    }
    list(out)
}

Orbit.combined <- function(t,state,parameters,verbose=FALSE){
        out <- rep(NA,6*Nt)
        ind <- 1:Nt
        xs <- state[(1:Nt)*6-5]
        ys <- state[(1:Nt)*6-4]
        zs <- state[(1:Nt)*6-3]
        for(j in 1:length(Ms)){
            out[(j-1)*6+1:3] <- state[(j-1)*6+4:6]
            R <- sqrt(sum(state[(j-1)*6+1:2]^2))
            z <- state[(j-1)*6+3]
            phi <- xy2phi(state[(j-1)*6+1],state[(j-1)*6+2])
            ar <- -numerical.derivative(function(R,z) phit(R,z,parameters),R,z)
            aphi <- 0
            az <- -numerical.derivative(function(z,R) phit(R,z,parameters),z,R)
            d3 <- ((xs[-j]-xs[j])^2+(ys[-j]-ys[j])^2+(zs[-j]-zs[j])^2)^1.5
            out[(j-1)*6+4] <- -mu*sum(Ms[-j]*(xs[j]-xs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6+ar*cos(phi)#kpc/Myr^2
            out[(j-1)*6+5] <- -mu*sum(Ms[-j]*(ys[j]-ys[-j])/d3)/kpc2au^2/kpcmyr2auyr*1e6+ar*sin(phi)
            out[(j-1)*6+6] <- -mu*sum(Ms[-j]*(zs[j]-zs[-j])/d3)/kpc2au^2/kpcmyr2auyr*1e6+az
        }
        list(out)
}

Orbit2.tide <- function(t,y,parameters){
    out <- rep(NA,6*Nt)
    ind <- 1:Nt
    xs <- y[(1:Nt)*6-5]
    ys <- y[(1:Nt)*6-4]
    zs <- y[(1:Nt)*6-3]
    for(j in 1:Nt){
        out[(j-1)*6+1:3] <- y[(j-1)*6+4:6]
        R <- sqrt(sum(y[(j-1)*6+1:2]^2))
        z <- y[(j-1)*6+3]
        phi <- xy2phi(y[(j-1)*6+1],y[(j-1)*6+2])
        if(gal.pot=='pmp'){
            ar <- pmp.R(R,z,parameters)
            az <- pmp.z(R,z,parameters)
        }else if(gal.pot=='hhmn'){
            ar <- hhmn.R(R,z,parameters)
            az <- hhmn.z(R,z,parameters)
        }
        aphi <- 0
        out[(j-1)*6+4] <- ar*cos(phi)#kpc/Myr^2
        out[(j-1)*6+5] <- ar*sin(phi)
        out[(j-1)*6+6] <- az
    }
    out
}
Orbit3.tide <- function(t,y,parameters){
    out <- rep(NA,6*Nt)
    ind <- 1:Nt
    xs <- y[(1:Nt)*6-5]
    ys <- y[(1:Nt)*6-4]
    zs <- y[(1:Nt)*6-3]
    for(j in 1:Nt){
        out[(j-1)*6+1:3] <- y[(j-1)*6+4:6]
        R <- sqrt(sum(y[(j-1)*6+1:2]^2))
        z <- y[(j-1)*6+3]
        phi <- xy2phi(y[(j-1)*6+1],y[(j-1)*6+2])
        if(gal.pot=='pmp'){
            ar <- pmp.R(R,z,parameters)
            az <- pmp.z(R,z,parameters)
        }else if(gal.pot=='hhmn'){
            ar <- hhmn.R(R,z,parameters)
            az <- hhmn.z(R,z,parameters)
        }
        aphi <- 0
        out[(j-1)*6+4] <- ar*cos(phi)#kpc/Myr^2
        out[(j-1)*6+5] <- ar*sin(phi)
        out[(j-1)*6+6] <- az
    }
    list(out)
}
Orbit.probability <- function(t,y,parameters){
###only the target star and the sun
    out <- y
    xs <- y[1]
    ys <- y[2]
    zs <- y[3]
    Ns <- length(xs)
    Rs <- sqrt(xs^2+ys^2)
    phis <- xy2phi(xs,ys)
    for(j in 1:Ns){
        out[(j-1)*6+1:3] <- y[(j-1)*6+4:6]
        R <- Rs[j]
        z <- zs[j]
        phi <- phis[j]
        if(gal.pot=='pmp'){
            ar <- pmp.R(R,z,parameters)
            az <- pmp.z(R,z,parameters)
        }else if(gal.pot=='hhmn'){
            ar <- hhmn.R(R,z,parameters)
            az <- hhmn.z(R,z,parameters)
        }
        aphi <- 0
        out[(j-1)*6+4] <- ar*cos(phi)#kpc/Myr^2
        out[(j-1)*6+5] <- ar*sin(phi)
        out[(j-1)*6+6] <- az
    }
####calculate new covariance matrix for Galacto/helio-centric orbit
    cov.old <- array(0,dim=c(6,6))
#cat('length(y)=',length(y),'\n')
    cov.old[upper.tri(cov.old,diag=TRUE)] <- y[7:length(y)]
    cov.old[lower.tri(cov.old)] <- t(cov.old)[lower.tri(cov.old)]
    dcov <- cal.dcov(Rs[1],phis[1],zs[1],parameters,cov.old)
    out[(Ns*6+1):length(out)] <- dcov[upper.tri(dcov,diag=TRUE)]
    out
}
cal.dcov <- function(R,phi,z,parameters,cov.old){
    dcov <- array(0,dim=dim(cov.old))
    JRz <- hhmn.calJRz(R,z,parameters)
    Jrr <- JRz2xyz(JRz,phi,R,z,parameters)
    ir <- 1:3
    iv <- 4:6
####covariance term
#    cat('diag(Cov[irt,iv])=',diag(Cov[irt,iv]),'\n')
#    cat('(Jrr)=',as.numeric(Jrr),'\n')
#    cat('(JRz)=',as.numeric(JRz),'\n')
    cov.ar <- Jrr%*%cov.old[ir,ir]
    cov.av <- Jrr%*%cov.old[ir,iv]#important for cov(r,r)
    dcov[ir,ir] <- 2*cov.old[ir,iv]
#    cat('dCov[ir,ir]=',dCov[ir,ir],'\n')
###cov(r,v)
    dcov[ir,iv] <- cov.old[iv,iv]+cov.ar
###cov(v,v)
    dcov[iv,iv] <- 2*cov.av
    return(dcov)
}
Orbit1.probability <- function(t,y,parameters){
###only the target star and the sun
    out <- y
    xs <- y[c(1,7)]
    ys <- y[c(2,8)]
    zs <- y[c(3,9)]
    xh <- diff(rev(xs))
    yh <- diff(rev(ys))
    zh <- diff(rev(zs))
    cov.old <- y[12:length(y)]
    Rs <- sqrt(xs^2+ys^2)
    phis <- xy2phi(xs,ys)
    Rh <- sqrt(xh^2+yh^2+zh^2)
    phih <- xy2phi(xh,yh)
    for(j in 1:2){
        out[(j-1)*6+1:3] <- y[(j-1)*6+4:6]
        R <- Rs[j]
        z <- zs[j]
        phi <- phis[j]
        if(gal.pot=='pmp'){
            ar <- pmp.R(R,z,parameters)
            az <- pmp.z(R,z,parameters)
        }else if(gal.pot=='hhmn'){
            ar <- hhmn.R(R,z,parameters)
            az <- hhmn.z(R,z,parameters)
        }
        aphi <- 0
        out[(j-1)*6+4] <- ar*cos(phi)#kpc/Myr^2
        out[(j-1)*6+5] <- ar*sin(phi)
        out[(j-1)*6+6] <- az
    }
####calculate new covariance matrix for Galacto/helio-centric orbit
    Npar <- sqrt(length(cov.old))
    cov.old <- t(array(cov.old,dim=c(Npar,Npar)))
    cov.old[lower.tri(cov.old)] <- cov.old[upper.tri(cov.old)]
    dcov <- array(0,dim(cov.old))
    ind1 <- (1:Npar)[-(7:12)]
    ind2 <- (1:Npar)[-(1:6)]
    cov.gal <- cov.old[ind1,ind1]
    cov.hel <- cov.old[ind2,ind2]
    dcov[ind1,ind1] <- calCov(Rs[1],phis[1],zs[1],cov.gal,parameters)
    dcov[ind2,ind2] <- calCov(Rh,phih,zh,cov.hel,parameters)
    out[-c(1:12)] <- as.numeric(dcov)
    list(out)
}

Orbit2.keplerian <- function(t,y){
    out <- rep(NA,6*Nt)
    mu <- 4*pi^2
    ind <- 1:Nt
    xs <- y[(1:Nt)*6-5]
    ys <- y[(1:Nt)*6-4]
    zs <- y[(1:Nt)*6-3]
    for(j in 1:Nt){
        out[(j-1)*6+1:3] <- y[(j-1)*6+4:6]
        d3 <- (xs[j]^2+ys[j]^2+zs[j]^2)^1.5
        a <- -mu*Mv*xs[j]/d3/kpc2au^2/kpcmyr2auyr*1e6
#        cat('a=',a,'\n')
#        cat('Mv=',Mv,'\n')
        out[(j-1)*6+4] <- a
        out[(j-1)*6+5] <- -mu*sum(Mv*ys[j]/d3)/kpc2au^2/kpcmyr2auyr*1e6
        out[(j-1)*6+6] <- -mu*sum(Mv*zs[j]/d3)/kpc2au^2/kpcmyr2auyr*1e6
    }
#    cat('head(out)=',head(out),'\n')
    out
}

Orbit2.combined <- function(t,y){
        out <- rep(NA,6*Nt)
        ind <- 1:Nt
        xs <- y[(1:Nt)*6-5]
        ys <- y[(1:Nt)*6-4]
        zs <- y[(1:Nt)*6-3]
        for(j in 1:length(Ms)){
            out[(j-1)*6+1:3] <- y[(j-1)*6+4:6]
            R <- sqrt(sum(y[(j-1)*6+1:2]^2))
            z <- y[(j-1)*6+3]
            phi <- xy2phi(y[(j-1)*6+1],y[(j-1)*6+2])
            ar <- -numerical.derivative(function(R,z) phit(R,z,parameters),R,z)
            aphi <- 0
            az <- -numerical.derivative(function(z,R) phit(R,z,parameters),z,R)
            d3 <- ((xs[-j]-xs[j])^2+(ys[-j]-ys[j])^2+(zs[-j]-zs[j])^2)^1.5
            out[(j-1)*6+4] <- -mu*sum(Ms[-j]*(xs[j]-xs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6+ar*cos(phi)#kpc/Myr^2
            out[(j-1)*6+5] <- -mu*sum(Ms[-j]*(ys[j]-ys[-j])/d3)/kpc2au^2/kpcmyr2auyr*1e6+ar*sin(phi)
            out[(j-1)*6+6] <- -mu*sum(Ms[-j]*(zs[j]-zs[-j])/d3)/kpc2au^2/kpcmyr2auyr*1e6+az
        }
        out
}

Orbit.separate <- function(t,state,parameters){
    nn <- length(state)/6
    out <- rep(NA,length(state))
    xs <- state[(1:(nn-1))*6-5]
    ys <- state[(1:(nn-1))*6-4]
    zs <- state[(1:(nn-1))*6-3]
    kpc2au <- 1e3*206265
    for(j in 1:length(xs)){
        out[(j-1)*6+1:3] <- state[(j-1)*6+4:6]
        d3 <- ((xs[-j]-xs[j])^2+(ys[-j]-ys[j])^2+(zs[-j]-zs[j])^2)^1.5
        R <- sqrt(sum(state[(j-1)*6+1:2]^2))
        z <- state[(j-1)*6+3]
        phi <- xy2phi(state[(j-1)*6+1],state[(j-1)*6+2])
        ar <- -numerical.derivative(function(R,z) phit(R,z,parameters),R,z)
        aphi <- 0
        az <- -numerical.derivative(function(z,R) phit(R,z,parameters),z,R)
                                        #            tmp <- mu*sum(Ms[-j]*(xs[j]-xs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6
                                        #            cat('tmp=',tmp,'\n')
                                        #            cat('ar=',ar,'\n')
        out[(j-1)*6+4] <- -mu*sum(Ms[-j]*(xs[j]-xs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6+ar*cos(phi)#kpc/Myr^2
        out[(j-1)*6+5] <- -mu*sum(Ms[-j]*(ys[j]-ys[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6+ar*sin(phi)
        out[(j-1)*6+6] <- -mu*sum(Ms[-j]*(zs[j]-zs[-j])/d3/kpc2au^2)/kpcmyr2auyr*1e6+az
    }
######calculate the state of companion separately
    x <- c(state[6*nn-5],state[6*nn-5]+(xs[1]-xs[2])*kpc2au)#au; separation from proxima and alpha Cen
    y <- c(state[6*nn-4],state[6*nn-4]+(ys[1]-ys[2])*kpc2au)
    z <- c(state[6*nn-3],state[6*nn-3]+(zs[1]-zs[2])*kpc2au)
    d3 <- (x^2+y^2+z^2)^1.5#au^3
#    cat('x=',x,'\n')
#    cat('y=',y,'\n')
#    cat('z=',z,'\n')
#    cat('d3=',d3,'\n')
    out[(nn-1)*6+1:3] <- state[(nn-1)*6+4:6]*1e6#au/yr; transfer yr to Myr time step
    out[(nn-1)*6+4] <- -mu*sum(Ms[2:1]*x/d3)*1e6#au/yr^2; transfer yr to Myr time step
    out[(nn-1)*6+5] <- -mu*sum(Ms[2:1]*y/d3)*1e6
    out[(nn-1)*6+6] <- -mu*sum(Ms[2:1]*z/d3)*1e6
    list(out)
}

Orbit.solar<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dR<- Rdot
    dRdot<- -numerical.derivative(function(R,z) phit(R,z,parameters),R,z)+R*phidot^2
    dphi<- phidot
    dphidot<- -2*Rdot/R*phidot
    dz<- zdot
    dzdot<- -numerical.derivative(function(z,R) phit(R,z,parameters),z,R)
    list(c(dR,dRdot,dphi,dphidot,dz,dzdot))
  })
}

Orbit.solar.gardner <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dR<- Rdot
    dRdot<- -numerical.derivative(function(R,z) Pott(R,z,parameters),R,z)+R*phidot^2
    dphi<- phidot
    dphidot<- -2*Rdot/R*phidot
    dz<- zdot
    dzdot<- -numerical.derivative(function(z,R) Pott(R,z,parameters),z,R)
    list(c(dR,dRdot,dphi,dphidot,dz,dzdot))
  })
}

tidal.force <- function(R,z,par.gal,par.tide){
  Vc <- Vcirc(R,par.gal)#kpc/myr
  dv.dr <- deriv.array(Vc,R)
  Omegac <- Vc/R
  A <- 1/2*(Omegac-dv.dr)
  B <- -1/2*(Omegac+dv.dr)
  A <- A/10^6 #Myr -> yr
  B <- B/10^6
  g1 <- -(A-B)*(3*A + B)
  g2 <- (A - B)^2
  rhot <- dent(R,z,par.gal)
  rhot.au <- rhot/10^9/206265^3
  g3 <- 4*pi*par.tide$mu*rhot.au-2*(B^2-A^2)
  return(list(G1=g1,G2=g2,G3=g3,A=A,B=B,rho=rhot))
}

xy2phi <- function(x,y){
    phi <- atan(y/x)
    inds <- which(x<0)
    phi[inds] <- phi[inds]+pi
    inds <- which(x>=0 & y<0)
    phi[inds] <- phi[inds]+2*pi
    return(phi)
}

xy2Rphi <- function(x,y){
    phi <- atan(y/x)
    inds <- which(x<0)
    phi[inds] <- phi[inds]+pi
    inds <- which(x>=0 & y<0)
    phi[inds] <- phi[inds]+2*pi
    R <- sqrt(x^2+y^2)
    return(cbind(R,phi))
}

###Note that the x,y,z here should be in the heliocentric frame, and x axis point to the GC while y axis point to the rotation direction
xyz2bl <- function(x,y,z){
  b <- rep(NA,length(x))
  l <- rep(NA,length(x))
  for(j in 1:length(x)){
    l.test <- atan(y[j]/x[j])
    if(x[j]<0){
      l[j] <- l.test+pi
    }else if((x[j]>0) & (y[j]<0)){
      l[j] <- l.test+2*pi
    }else{
      l[j] <- l.test
    }
    b[j] <- atan(z[j]/sqrt(x[j]^2+y[j]^2))
  }
  return(cbind(b,l))
}

###Note that the x,y,z here should be in the heliocentric frame, and x axis point to the GC while y axis point to the rotation direction
xyz2bl.vec <- function(x,y,z){
    b <- atan(z/sqrt(x^2+y^2))
    l <- atan(y/x)
    inds <- which(x<0)
    l[inds] <- l[inds]+pi
    inds <- which(x>=0 & y<0)
    l[inds] <- l[inds]+2*pi
    return(cbind(b,l))
}

xyz2equ <- function(x,y,z){
####cal. delta
    delta <- atan(z/sqrt(x^2+y^2))
###cal. alpha
    alpha <- atan(y/x)
    inds <- which(x<0)
    alpha[inds] <- alpha[inds]+pi
    inds <- which(x>=0 & y<0)
    alpha[inds] <- alpha[inds]+2*pi
    return(cbind(alpha,delta))
}

angle.between <- function(a,b){
  ra <- sqrt(a[1]^2+a[2]^2+a[3]^2)
  rb <- sqrt(b[1]^2+b[2]^2+b[3]^2)
  costh <- a%*%b/(ra*rb)
  sinth <- sqrt(1-cos(costh)^2)
  return(c(cos=costh,sin=sinth))
}

transform2apex <- function(b,l,b.sun,l.sun){
  b.apex <- rep(NA,length(b))
  l.apex <- rep(NA,length(b))
  x.apex <- rep(NA,length(b))
  y.apex <- rep(NA,length(b))
  z.apex <- rep(NA,length(b))
  x <- cos(b)*cos(l)
  y <- cos(b)*sin(l)
  z <- sin(b)
  x.sun <- cos(b.sun)*cos(l.sun)
  y.sun <- cos(b.sun)*sin(l.sun)
  z.sun <- sin(b.sun)
  r.sun <- c(x.sun,y.sun,z.sun)
  for(j in 1:length(b)){
    rotz <- matrix(data=c(cos(l.sun),sin(l.sun),0.0,-sin(l.sun),cos(l.sun),0.0,0.0,0.0,1.0),nrow=3,ncol=3,byrow=TRUE)#o-xyz -> o-x'y'z'
    #roty <- matrix(data=c(cos(b.sun),0.0,sin(b.sun),0.0,1.0,0.0,-sin(b.sun),0.0,cos(b.sun)),nrow=3,ncol=3,byrow=TRUE)#o-x'y'z' ->o-x''y''z''
    roty <- matrix(data=c(sin(b.sun),0.0,-cos(b.sun),0.0,1.0,0.0,cos(b.sun),0.0,sin(b.sun)),nrow=3,ncol=3,byrow=TRUE)#o-x'y'z' ->o-x''y''z''; NGP -> apex
    r.apex <- roty%*%rotz%*%c(x[j],y[j],z[j])
    x.apex[j] <- r.apex[1]
    y.apex[j] <- r.apex[2]
    z.apex[j] <- r.apex[3]
    l.test <- atan(y.apex[j]/x.apex[j])
    if(x.apex[j]<0){
      l.apex[j] <- l.test+pi
    }else if((x.apex[j]>0) & (y.apex[j]<0)){
      l.apex[j] <- l.test+2*pi
    }else{
      l.apex[j] <- l.test
    }
    b.apex[j] <- atan(z.apex[j]/sqrt(x.apex[j]^2+y.apex[j]^2))
  }
  return(cbind(b.apex,l.apex,x.apex,y.apex,z.apex))
}

DirGen <- function(x.sun,y.sun,z.sun){
#vertical vector in xy plane
  r.sun <- c(x.sun,y.sun,z.sun)
  la <- atan(-x.sun/y.sun)
  ba <- 0
  xa <- cos(la)
  ya <- sin(la)
  za <- 0
  ra <- c(xa,ya,za)
#vector which is perpendicular with the above two vectors
  rd <- cross(ra,r.sun)
  xd <- rd[1]
  yd <- rd[2]
  zd <- rd[3]
###############generate vectors which are perpendicular to solar apex
  theta <- seq(0,2*pi,by=0.01)
  d <- array(data=NA,dim=c(length(theta),3))
  for(i in 1:length(theta)){
    d[i,] <- cos(theta[i])*rd+sin(theta[i])*ra
  }
  bl.vert <- xyz2bl.vec(d[,1],d[,2],d[,3])
  return(bl.vert)
}

invariant2galactic <- function(x.inv,y.inv,z.inv,b.inv,l.inv){
  x <- rep(NA,length(x.inv))
  y <- rep(NA,length(x.inv))
  z <- rep(NA,length(x.inv))
#  cat("x':",x.inv,"y':",y.inv,"z':",z.inv,"\n")
  for(j in 1:length(x.inv)){
    rotz <- matrix(data=c(cos(l.inv),-sin(l.inv),0.0,sin(l.inv),cos(l.inv),0.0,0.0,0.0,1.0),nrow=3,ncol=3,byrow=TRUE)#o-xyz -> o-x'y'z'
    roty <- matrix(data=c(sin(b.inv),0.0,cos(b.inv),0.0,1.0,0.0,-cos(b.inv),0.0,sin(b.inv)),nrow=3,ncol=3,byrow=TRUE)#o-x'y'z' ->o-x''y''z''; NGP -> NEP
    #roty <- matrix(data=c(sin(b.inv),0.0,-cos(b.inv),0.0,1.0,0.0,cos(b.inv),0.0,sin(b.inv)),nrow=3,ncol=3,byrow=TRUE)#o-x'y'z' ->o-x''y''z''
    r <- rotz%*%roty%*%c(x.inv,y.inv,z.inv)
    x[j] <- r[1]
    y[j] <- r[2]
    z[j] <- r[3]
  }
  return(cbind(x,y,z))
}

#########transfer from a heliocentric coordinate system to a galactic centric coordinate
hel2gal<- function(Dc, b, l, Rs,Zs){
    Xh <- -Dc*cos(b)*cos(l)
    Yh <- -Dc*cos(b)*sin(l)
    Zh <- Dc*sin(b)
    Xg <- Rs+Xh
    Yg <- Yh
    Rg <- sqrt(Xg^2+Yg^2)
    Zg <- Zs+Zh
    phig<- atan(Yg/Xg)
    for(i in 1:length(phig)){
        if(Xg[i]>=0){
            phig[i]<- phig[i]
        }else if(Xg[i]<0 & Yg[i]<=0){
            phig[i] <- phig[i]-pi
        }else{
            phig[i] <- phig[i]+pi
        }
    }
    return(list(phi=phig,R=Rg,X=Xg,Y=Yg, Z=Zg,Xh=Xh,Yh=Yh,Zh=Zh))
}
###EL calculation
calEL <- function(out){
    Mv <- 1.32e12#Galaxy as a point-like object; Keplerian potential
    r <- sqrt(rowSums(out[,1:3,drop=FALSE]^2))#kpc
    v <- sqrt(rowSums(out[,4:6,drop=FALSE]^2))#kpc/myr
    K <- 0.5*(v*kpcmyr2auyr)^2
    V <- -(2*pi)^2*Mv/(r*kpc2au)
    E <- K+V
    dE <- (E-E[1])/E[1]
    ##angular momentum
    x <- out[,1]
    y <- out[,2]
    z <- out[,3]
    vx <- out[,4]
    vy <- out[,5]
    vz <- out[,6]
    L <- sqrt((x*vy-y*vx)^2+(y*vz-z*vy)^2+(z*vx-x*vz)^2)
    Lz <- y*vx-x*vy
    dLz <- (Lz-Lz[1])/Lz[1]
    dL <- (L-L[1])/L[1]
    return(list(E=E,dE=dE,L=L,dL=dL,Lz=Lz,dLz=dLz))
}
subsampling <- function(var,N=10){
    if(is.null(dim(var))){
        vmid <- seq(min(var),max(var),length.out=N+1)
        inds <- sapply(1:N, function(j) which.min(abs(var-vmid[j])))
    }else{
        vars <- list()
        for(j in 1:ncol(var)){
            vars[[j]] <- seq(min(var[,j]),max(var[,j]),length.out=N)
        }
        inds <- c()
        grid <- expand.grid(vars)
        for(j in 1:nrow(grid)){
            tmp <- t(t(var)-as.numeric(grid[j,]))
            inds <- c(inds,which.min(rowSums(tmp^2)))
        }
    }
    return(inds)
}

subsampling3D <- function(var,N=10){
    var1 <- seq(min(var[,1]),max(var[,1]),length.out=N+1)
    var2 <- seq(min(var[,2]),max(var[,2]),length.out=N+1)
    D3 <- FALSE
    if(ncol(var)==3){
        D3 <- TRUE
    }
    if(D3){
        var3 <- seq(min(var[,3]),max(var[,3]),length.out=N+1)
    }
    inds <- c()
    for(i in 1:N){
        for(j in 1:N){
            if(D3){
                N3 <- N
            }else{
                N3 <- 1
            }
            for(k in 1:N3){
                l1 <- which(var[,1]>=var1[i] & var[,1]<=var1[i+1])
                l2 <- which(var[,2]>=var2[j] & var[,2]<=var2[j+1])
                ind12 <- intersect(l1,l2)
                if(D3){
                    l3 <- which(var[,3]>=var3[k] & var[,3]<=var3[k+1])
                    ind123 <- intersect(intersect(l1,l2),l3)
                }
                if(D3){
                    if(length(ind123)>0){
                        ind <- sample(ind123,1)
                    }else if(length(ind12)>0){
                        ind <- sample(ind12,1)
                    }else{
                        ind <- sample(l1,1)
                    }
                }else{
                    if(length(ind12)>0){
                        ind <- sample(ind12,1)
                    }else if(length(l1)>0){
                        ind <- sample(l1,1)
                    }else{
                        ind <- c()
                    }
                }
                inds <- c(inds,ind)
            }
        }
    }
    return(inds)
}

#for 1D density
subsampling1 <- function(var,N=10){
    if(is.null(dim(var))){
        fit <- density(var)
        fit.fun <- approxfun(fit$x,fit$y)
        inds <- sort(sample(length(var), N, prob=1/fit.fun(var)))
    }
    return(inds)
}
#use kde for high dimension kernel density
subsampling2 <- function(var,N=1e3){
    ker <- kde(data,eval.points=data)
    inds <- sample(1:nrow(data), 10*N, prob=1/ker$estimate)
    inds <- sort(unique(inds))[1:N]
    return(inds)
}
##time consuming
subsampling3D.cycle <- function(var,N=1e3){
    fit1 <- density(var[,1])
    fit1.fun <- approxfun(fit1$x,fit1$y)
    fit2 <- density(var[,2])
    fit2.fun <- approxfun(fit2$x,fit2$y)
    D3 <- FALSE
    if(ncol(var)==3) D3 <- TRUE
    if(D3){
        fit3 <- density(var[,3])
        fit3.fun <- approxfun(fit3$x,fit3$y)
    }
    Ntry <- 1e3
    for(j in 1:Ntry){
        ind1 <- sample(nrow(var), N, prob=1/fit1.fun(var[,1]))
        ind2 <- sample(nrow(var), N, prob=1/fit2.fun(var[,2]))
        if(D3){
            ind3 <- sample(nrow(var), N, prob=1/fit3.fun(var[,3]))
            ind <- intersect(intersect(ind1,ind2),ind3)
        }else{
            ind <- intersect(ind1,ind2)
        }
        if(length(ind)>0){
            inds <- unique(c(inds,ind))
        }
        cat('length(inds)=',length(inds),'\n')
        if(length(inds)>N) break()
    }
    return(inds[1:N])
}
####transform velocity in equatorial coordinates to velocity in galactic coordinates
####units: pm.ra and pm.de: mas/yr; rv: km/s; ra, de: rad; distance from the target to the Sun d: kpc
####UVW:helio-static Galactic velocity; U(positive toward the Galactic Center), V(postive toward the Galactic rotation), W(positive toward to the North Galactic pole)
Ve2hg <- function(pm.ra,pm.de,rv,ra,de,d,test=FALSE){
    auyr2kms <- 4.74047
    vra <- pm.ra*d#au/yr
    vde <- pm.de*d#au/yr
    vp <- sqrt(vra^2+vde^2)
    vr <- rv/auyr2kms#au/yr
    vx.equ <- vr*cos(de)*cos(ra)-vde*sin(de)*cos(ra)-vra*sin(ra)##note: vr is positive if the star is moving away from the Sun
    vy.equ <- vr*cos(de)*sin(ra)-vde*sin(de)*sin(ra)+vra*cos(ra)
    vz.equ <- vr*sin(de)+vde*cos(de)
    v.equ <- cbind(vx.equ,vy.equ,vz.equ)
    vhg <- e2g.vel(vx.equ,vy.equ,vz.equ)
    if(test){
        bl <- equ2gal(ra[1],de[1])
        ve <- vhg[1,]
        re <- bl2xyz(bl[1],bl[2])
        cat('Galactic distance: ',d[1]*re,'pc\n')
        cat('Galactic velocity: ',vhg[1,],'au/yr\n')
        cat('r1*v1=',d[1]*sum(re*vhg[1,]),'\n')
        vre <- sum(ve*re)
        vpe <- sqrt(ve%*%ve-vre^2)
        cat('vp=',vp[1],'; vpe=',vpe,'\n')
        cat('vr=',vr[1],'; vre=',vre,'\n')
    }
    return(list(ve=v.equ,vhg=vhg,vp=vp))#in unit of au/yr
}

####transform helio-static velocity in the Galactic coordinates into proper motions
###ra and dec are in unit of rad
hg2pm <- function(vx,vy,vz,alpha,delta,d){
    ve <- g2e.vel(vx,vy,vz)
    vra <- ve[,2]*cos(alpha)-ve[,1]*sin(alpha)
    vxy <- ve[,1]*cos(alpha)+ve[,2]*sin(alpha)
    vde <- ve[,3]*cos(delta)-vxy*sin(delta)
    vr <- ve[,3]*sin(delta)+vxy*cos(delta)#au/yr
    pmra <- vra/d#mas/yr
    pmde <- vde/d
    return(cbind(vra,vde,vr,pmra,pmde))
}

####transform from static Galactic reference frame to rotation Galactic reference frame, or from vx, vy, vz to U, V, W
uvw2gal <- function(vx,vy,vz,phi){
    U <- -(vx*cos(phi)+vy*sin(phi))
    V <- vx*sin(phi)-vy*cos(phi)
    W <- vz
    return(cbind(U,V,W))
}


#####convertions according to http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1987AJ.....93..864J&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
###


#version final: adapt integrator.R to integrate TNO orbit
tide.HM<-function(t,state,par.tide){
with(as.list(c(state,par.tide)),{
  Omegar <- Omega-Omega0*t
  dL <- 0
  dG <- -5*L^2*(L^2-G^2)/(2*mu^2)*(cos(w)*sin(w)*(g3*(1-Theta^2/G^2)+(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*Theta^2/G^2-g1*cos(Omegar)^2-g2*sin(Omegar)^2)-(g1-g2)*(cos(w)^2-sin(w)^2)*cos(Omegar)*sin(Omegar)*Theta/G)
  dTheta <- L^2*(g1-g2)/(2*mu^2)*(5*(L^2-G^2)*Theta/G*cos(w)*sin(w)*(cos(Omegar)^2-sin(Omegar)^2)+sin(Omegar)*cos(Omegar)*(G^2-Theta^2+5*(L^2-G^2)*(cos(w)^2-sin(w)^2*Theta^2/G^2)))
  dw <- L^2*G/(2*mu^2)*(g3*(1-5*sin(w)^2*(1-L^2*Theta^2/G^4))+(g1*cos(Omegar)^2+g2*sin(Omegar)^2)*(1-5*cos(w)^2)-5*(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*L^2*Theta^2/G^4*sin(w)^2+5*(g1-g2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar)*(G^2+L^2)*Theta/G^3)
  dOmega <- L^2/(2*G*mu^2)*((g1*sin(Omegar)^2+g2*cos(Omegar)^2-g3)*(G^2+5*(L^2-G^2)*sin(w)^2)*Theta/G-5*(g1-g2)*(L^2-G^2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar))
  return(list(c(dw,dOmega,dL,dG,dTheta)))
})
}

tide.CM<-function(t,state,par.tide){
with(as.list(c(state,par.tide)),{
  mu <- par.tide$mu
  g1 <- g1.func(t)
  g2 <- g2.func(t)
  g3 <- g3.func(t)
#  g1 <- par.tide$g1
#  g2 <- par.tide$g2
#  g3 <- par.tide$g3
  Omega0 <- par.tide$Omega0
  r <- sqrt(x^2+y^2+z^2)
  xp <- x*cos(Omega0*t)+y*sin(Omega0*t)
  yp <- -x*sin(Omega0*t)+y*cos(Omega0*t)
  dx <- xdot
  dy <- ydot
  dz <- zdot
  dxdot <- -mu*x/r^3-g1*xp*cos(Omega0*t)+g2*yp*sin(Omega0*t)
  dydot <- -mu*y/r^3-g1*xp*sin(Omega0*t)-g2*yp*cos(Omega0*t)
  dzdot <- -mu*z/r^3-g3*z
  return(list(c(dx,dy,dz,dxdot,dydot,dzdot)))
})
}

map3 <- function(t,dt,state,par.tide){
  L <- state[1]
  G <- state[2]
  Theta <- state[3]
  w <- state[4]
  Omega <- state[5]
  Omega0 <- par.tide$Omega0
  mu <- par.tide$mu
  g1 <- par.tide$g1
  g2 <- par.tide$g2
  g3 <- par.tide$g3
#  g1 <- par.tide$g1
#  g2 <- par.tide$g2
#  g3 <- par.tide$g3
  Omegar <- Omega-Omega0*t
  dL <- 0
  dG <- -5*L^2*(L^2-G^2)/(2*mu^2)*(cos(w)*sin(w)*(g3*(1-Theta^2/G^2)+(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*Theta^2/G^2-g1*cos(Omegar)^2-g2*sin(Omegar)^2)-(g1-g2)*(cos(w)^2-sin(w)^2)*cos(Omegar)*sin(Omegar)*Theta/G)
  dTheta <- L^2*(g1-g2)/(2*mu^2)*(5*(L^2-G^2)*Theta/G*cos(w)*sin(w)*(cos(Omegar)^2-sin(Omegar)^2)+sin(Omegar)*cos(Omegar)*(G^2-Theta^2+5*(L^2-G^2)*(cos(w)^2-sin(w)^2*Theta^2/G^2)))
  dw <- L^2*G/(2*mu^2)*(g3*(1-5*sin(w)^2*(1-L^2*Theta^2/G^4))+(g1*cos(Omegar)^2+g2*sin(Omegar)^2)*(1-5*cos(w)^2)-5*(g1*sin(Omegar)^2+g2*cos(Omegar)^2)*L^2*Theta^2/G^4*sin(w)^2+5*(g1-g2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar)*(G^2+L^2)*Theta/G^3)
  dOmega <- L^2/(2*G*mu^2)*((g1*sin(Omegar)^2+g2*cos(Omegar)^2-g3)*(G^2+5*(L^2-G^2)*sin(w)^2)*Theta/G-5*(g1-g2)*(L^2-G^2)*cos(w)*sin(w)*cos(Omegar)*sin(Omegar))
  dOmegar <- dOmega-Omega0
######calculate the second derivatives
  x <- c(L,G,Theta,w,Omegar)
  f1 <- function(x) c(dL, dG, dTheta, dw, dOmegar)
  f2 <- function(x) jacobian(f1,c(L,G,Theta,w, Omegar))%*%f1(x)
######calculate the third derivatives
  f3 <- function(x) jacobian(f2,c(L,G,Theta,w, Omegar))%*%f1(x)
  L <- L+dt*f1(x)[1]+dt^2/2*f2(x)[1]+dt^3/6*f3(x)[1]
  G <- G+dt*f1(x)[2]+dt^2/2*f2(x)[2]+dt^3/6*f3(x)[2]
  Theta <- Theta+dt*f1(x)[3]+dt^2/2*f2(x)[3]+dt^3/6*f3(x)[3]
  w <- w+dt*f1(x)[4]+dt^2/2*f2(x)[4]+dt^3/6*f3(x)[4]
  Omega <- Omega+dt*f1(x)[5]+dt^2/2*f2(x)[5]+dt^3/6*f3(x)[5]
  return(c(t+dt,w,Omega,L,G,Theta))
}

############################################################
##The following functions are for LARKS model of Oort cloud
############################################################
cv2ks <- function(x,y,z,vx,vy,vz,alpha){
  r <- sqrt(x^2+y^2+z^2)
  if(x>=0){
    u <- sqrt(alpha/(2*(r+x)))*c(0,r+x,y,z)
  }else{
    u <- sqrt(alpha/(2*(r-x)))*c(-z,y,r-x,0)}
  U <- 2/alpha*c(u[1]*vx+u[4]*vy-u[3]*vz,u[2]*vx+u[3]*vy+u[4]*vz,-u[3]*vx+u[2]*vy-u[1]*vz,-u[4]*vx+u[1]*vy+u[2]*vz)
  val <- u[2]*U[1]-u[1]*U[2]-u[4]*U[3]+u[3]*U[4]
###give a test of this identity
  if(abs(val)>0.001){
     print("u1*U0-u0*U1-u3*U2+u2*U3!=0")}
  return(c(u=u,U=U))
}

ks2cv <- function(u,U,alpha){
  x <- (u[1]^2+u[2]^2-u[3]^2-u[4]^2)/alpha
  y <- 2*(u[2]*u[3]+u[1]*u[4])/alpha
  z <- 2*(u[2]*u[4]-u[1]*u[3])/alpha
  r <- sqrt(x^2+y^2+z^2)
  vx <- 1/(2*r)*(u[1]*U[1]+u[2]*U[2]-u[3]*U[3]-u[4]*U[4])
  vy <- 1/(2*r)*(u[4]*U[1]+u[3]*U[2]+u[2]*U[3]+u[1]*U[4])
  vz <- 1/(2*r)*(-u[3]*U[1]+u[4]*U[2]-u[1]*U[3]+u[2]*U[4])
  return(c(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz))
}

H0 <- function(x,y,z,vx,vy,vz,par.tide){
  mu <- par.tide$mu
  r <- sqrt(x^2+y^2+z^2)
  val <- 1/2*(vx^2+vy^2+vz^2)-mu/r
  return(val)
}

H1 <- function(t,x,y,z,par.tide){
  g1 <- par.tide$g1
  g2 <- par.tide$g2
  g3 <- par.tide$g3
  Omega0 <- par.tide$Omega0
  xp <- x*cos(Omega0*t)+y*sin(Omega0*t)
  yp <- -x*sin(Omega0*t) + y*cos(Omega0*t)
  val <- 1/2*(g1*xp^2+g2*yp^2+g3*z^2)
  return(val)
}

pH1.pu <- function(t,x,y,z,par.tide){
  r <- sqrt(x^2+y^2+z^2)
  if(x>=0){
    u <- sqrt(alpha/(2*(r+x)))*c(0,r+x,y,z)
  }else{
    u <- sqrt(alpha/(2*(r-x)))*c(-z,y,r-x,0)}
  val <- grad(function(u) H1(t,x,y,z,par.tide),u)
}

K0 <- function(u,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)[1]
  y <- ks2cv(u,U,alpha)[2]
  z <- ks2cv(u,U,alpha)[3]
  vx <- ks2cv(u,U,alpha)[4]
  vy <- ks2cv(u,U,alpha)[5]
  vz <- ks2cv(u,U,alpha)[6]
  kappa0 <- H0(x,y,z,vx,vy,vz,par.tide)
  return(kappa0)
}

K1 <- function(u,t,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)["x"]
  y <- ks2cv(u,U,alpha)["y"]
  z <- ks2cv(u,U,alpha)["z"]
  kappa1 <- H1(t,x,y,z,par.tide)
  return(kappa1)
}

Us.ks <- function(t,u,U,alpha,par.tide){
  -K0(u,U,alpha,par.tide)-K1(u,t,U,alpha,par.tide)
}

Us.cv <- function(t,x,y,z,vx,vy,vz,par.tide){
  -H0(x,y,z,vx,vy,vz,par.tide)-H1(t,x,y,z,par.tide)
}

M0 <- function(u,U,alpha,par.tide){
  1/2*sum(U^2)+(4*Us/alpha^2)*sum(u^2)
}

M1 <- function(t,x,y,z,par.tide){
  4*sum(u^2)/alpha^2*H1(t,x,y,z,par.tide)
}

M1.ts <- function(t,u,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)[1]
  y <- ks2cv(u,U,alpha)[2]
  z <- ks2cv(u,U,alpha)[3]
  val <- M1(t,x,y,z,par.tide)
  return(val)
}

hess.M1 <- function(t,u,U,alpha,par.tide){
  x <- ks2cv(u,U,alpha)[1]
  y <- ks2cv(u,U,alpha)[2]
  z <- ks2cv(u,U,alpha)[3]
  val <- hessian(function(u) M1.ts(t,u,U,par.tide),u)
  return(val)
}

Phi0 <- function(u,t,U,Us,alpha,Delta,par.tide){
  if(Us>0){
    w <- 2*sqrt(2*Us)/alpha
    v <- u*cos(w*Delta)+U*w^-1*sin(w*Delta)
    V <- -u*w*sin(w*Delta)+U*cos(w*Delta)
    Vs <- Us
    t.new <- t+2*Delta/alpha^2*(sum(u^2)+sum(U^2)/w^2)+2*(u%*%U-v%*%V)/(alpha^2*w^2)
    C1 <- u%*%u + U%*%U/w^2
    C2 <- v%*%v + V%*%V/w^2
    err <- abs((C1-C2)/C1)
    if(err>0.01){
      print("The expression 'u^2+U^2*w^-2' is not invariant!")
      print(err)}
  }else{
      w <- 2*sqrt(-2*Us)/alpha
      v <- u*cosh(w*Delta)+U*w^-1*sinh(w*Delta)
      V <- u*w*sinh(w*Delta)+U*cosh(w*Delta)
      Vs <- Us
      t.new <- t+2*Delta/alpha^2*(sum(u^2)-sum(U^2)/w^2)-2*(u%*%U-v%*%V)/(alpha^2*w^2)
    C1 <- u%*%u - U%*%U/w^2
    C2 <- v%*%v - V%*%V/w^2
    err <- abs((C1-C2)/C1)
    if(err>0.01){
      print("The expression 'u^2+U^2*w^-2' is not invariant!")
      print(err)}
    }
  #GIVE A TEST FOR U AND u
  return(c(u=v,t=t.new,U=V,Us=Vs))
}

F <- function(u,t,U,alpha,par.tide){
  8*K1(u,t,U,alpha,par.tide)/alpha^2*u + 4*sum(u^2)/alpha^2*grad(function(u) K1(u,t,U,alpha,par.tide),u)
}

Fs <- function(u,t,U,alpha,par.tide){
  4*sum(u^2)/alpha^2*grad(function(t) K1(u,t,U,alpha,par.tide),t)
}

Phi1 <- function(u,t,U,Us,alpha,Delta,par.tide){
  v <- u
  t.new <- t
  V <- U-Delta*F(u,t,U,alpha,par.tide)
  Vs <- Us-Delta*Fs(u,t,U,alpha,par.tide)
  return(c(u=v,t=t.new,U=V,Us=Vs))
}

Phic <- function(u,t,U,Us,alpha,Delta,par.tide){
  v <- u
  t.new <- t
  V <- U-2*Delta*jacobian(function(u) F(u,t,U,alpha,par.tide), u)%*%F(u,t,U,alpha,par.tide)
  Vs <- Us-2*Delta*sum(jacobian(function(t) F(u,t,U,alpha,par.tide), t)*F(u,t,U,alpha,par.tide))
  return(c(u=v,t=t.new,U=V,Us=Vs))
}

Phih <- function(u,t,U,Us,alpha,h,par.tide){
  d1 <- h/12
  d2 <- 5/12*h
  c2 <- (1/2-sqrt(5)/10)*h
  c3 <- h/sqrt(5)
  q <- -h^3*(13-5*sqrt(5))/288
  D <- c(q,d1,c2,d2,c3,d2,c2,d1,q)
  funs <- list(Phic,Phi1,Phi0,Phi1,Phi0,Phi1,Phi0,Phi1,Phic)
  for(k in 1:3){
    Delta <- D[k]
    f <- funs[[k]]
    map <- as.numeric(f(u,t,U,Us,alpha,Delta,par.tide))
    u <- map[1:4]
    t <- map[5]
    U <- map[6:9]
    Us <- map[10]
  }
  return(c(u=u,t=t,U=U,Us=Us))
}

#Cartesian to Kepler; ref. http://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
cs2oe <- function(cs,Mc=1){
  R <- cs[1:3]
  V <- cs[4:6]
  r <- sqrt(R%*%R)
  mu <- Mc*4*pi^2#6.67384e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2
  v <- sqrt(V%*%V)
######equinoctial orbital element
  H <- cross(R,V)
  h <- sqrt(sum(H^2))
#  cat('dim(V)=',dim(V),'; dim(H)=',dim(H),'\n')
  E <- cross(V,H)/mu-R/r
  e <- sqrt(sum(E^2))#1 kepler: eccentricity
  p <- sum(h^2)/mu
  a <- h^2/(mu*(1-e^2))#2 kepler: semi-major axis
  i <- acos(H[3]/h)#3 kepler: inclination 0<=i<pi
  N <- cross(c(0,0,1),H)
  n <- sqrt(sum(N^2))
  if(N[2]<0){
    Omega <- 2*pi-acos(N[1]/n)
  }else{
    Omega <- acos(N[1]/n)
  }  #4 kepler: Right ascensioin of ascending node
  if(E[3]<0){
    w <- 2*pi-acos(sum(N*E)/(n*e))
  }else{
    w <- acos(sum(N*E)/(n*e))
  }#5 kepler: true anomaly
  p <- a*(1-e^2)
  nu <- atan2(sqrt(p/mu)*R%*%V,p-r)
  Ea <- 2*atan(tan(nu/2)*sqrt((1+e)/(1-e)))
  M <- Ea-e*sin(Ea)#6 kepler: mean anomaly
  q <- a*(1-e)
  oe <- c(e,q,i,w,Omega,M)
  return(oe)
}

###ref. https://downloads.rene-schwarz.com/download/M002-Cartesian_State_Vectors_to_Keplerian_Orbit_Elements.pdf
#and http://ccar.colorado.edu/asen5070/handouts/cart2kep2002.pdf
cs2oe.vec <- function(cs,Mc=1,t=NULL){
#  mu <- Mc*6.67428e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2
  mu <- Mc*4*pi^2
  if(is.null(dim(cs))) cs <- matrix(cs,nrow=1)
  Norb <- nrow(cs)
  R <- cs[,1:3,drop=FALSE]
  V <- cs[,4:6,drop=FALSE]
  r <- sqrt(rowSums(R*R))
  v <- sqrt(rowSums(V*V))
######equinoctial orbital element
  H <- cross(R,V)#angular momentum
  if(is.null(dim(H))){
      H <- matrix(H,nrow=1)
  }
  h <- sqrt(rowSums(H^2))#
###
  evec <- cross(V,H)/mu-R/r#eccentricity vector
  nvec <- cross(matrix(rep(c(0,0,1),Norb),byrow=TRUE,nrow=Norb),H)
  if(is.null(dim(evec))){
      evec <- matrix(evec,nrow=1)
      nvec <- matrix(nvec,nrow=1)
  }
  e <- sqrt(rowSums(evec^2))
##semi-major axis
  a <- 1/(2/r-v^2/mu)
###periapsis
  q <- a*(1-e)
###true anomaly
  p <- a*(1-e^2)
  rv <- rowSums(R*V)
  nu <- atan2(sqrt(p/mu)*rv,p-r)
###eccentric anomaly
  E <- 2*atan(tan(nu/2)*sqrt((1-e)/(1+e)))
###inclination
  i <- acos(H[,3]/h)#3 kepler: inclination 0<=i<pi
##longitude of the ascending node; Omega
  Omega <- atan2(H[,1],-H[,2])
  n <- sqrt(rowSums(nvec^2))
##argument of latitude
  w <- rep(0,nrow(nvec))
  ne <- rowSums(nvec*evec)
  inds <- which(n!=0&e!=0)
  w[inds] <- acos(ne[inds]/(n[inds]*e[inds]))
  inds <- which(evec[,3]<0)
  w[inds] <- 2*pi-w[inds]
##mean anomaly M
  M <- E-e*sin(E)
  oe <- cbind(e,q,i,w,Omega,M,a,nu,E)
  return(oe)
}

eps3 <- function(m,e,x){
    t1 <- cos(x)
    t2 <- -1+e*t1
    t3 <- sin(x)
    t4 <- e*t3
    t5 <- -x+t4+m
    t6 <- t5/(1/2*t5*t4/t2+t2)
    return(t5/((1/2*t3-1/6*t1*t6)*e*t6+t2))
}
KeplerStart3 <- function(m,e){
    t34 <- e^2
    t35 <- e*t34
    t33 <- cos(m)
    return(m+(-1/2*t35+e+(t34+3/2*t33*t35)*t33)*sin(m))
}
kep.murison2 <- function(m,e,tol=1e-6){
    Mnorm <- m%%(2*pi)
    E0 <- KeplerStart3(Mnorm,e)
    Ntry <- 1000
    for(k in 1:Ntry){
        E1 <- E0-eps3(Mnorm,e,E0)
        if(k==Ntry) 'Kepler solver failed to converge!\n'
        if(abs(E1-E0)<tol) break()
        E0 <- E1
    }
    return(E1)
}

kep.nr <- function(M,e,tol=1e-6){
    E <- E0 <- M
    Ntry <- 1e6
    for(j in 1:Ntry){
        E <- E0-(E0-e*sin(E0)-M)/(1-e*cos(E0))
        if(abs(E-E0)<tol) break()
    }
    return(E)
}
kep.mt <- function(m,e,tol=1e-6){
    E <- rep(NA,length(m))
    Ntt <- 10000
    for(j in 1:length(m)){
        E0 <- m[j]
        for(k in 1:Ntt){
            E1 = E0-(E0-e[j]*sin(E0)-m[j])/(sqrt((1-e[j]*cos(E0))^2-(E0-e[j]*sin(E0)-m[j])*(e[j]*sin(E0))))
            if(abs(E1-E0)<tol) break()
            if(j==Ntt) cat('Keplerian solver does not converge!\n')
            E0 <- E1
        }
        E[j] <- E1%%(2*pi)
    }
    return(E)
}

###https://downloads.rene-schwarz.com/download/M001-Keplerian_Orbit_Elements_to_Cartesian_State_Vectors.pdf
oe2cs.vec <- function(oe,Mc=1){
  #deriv Ea by solving Kepler's equation M=E-e*sinE for E
  e <- oe[,1]
  q <- oe[,2]
  i <- oe[,3]
  w <- oe[,4]
  Omega <- oe[,5]
  M <- oe[,6]
  a <- q/(1-e)
  mu <- 6.67428e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2#m^3/s^2
  mu <- Mc*mu#
#  Ea <- kep.nr(M,e,tol=1e-10)#vectorized
  Ea <- kep.mt(M,e,tol=1e-10)
#  Ea <- kep.murison2(M,e,tol=1e-10)

##true anomaly
  nu <- 2*atan2(sqrt(1+e)*sin(Ea/2),sqrt(1-e)*cos(Ea/2))
  r <- a*(1-e*cos(Ea))
  o <- cbind(r*cos(nu),r*sin(nu),0)
  v <- sqrt(mu*a)/r
  vo <- cbind(-v*sin(Ea),v*sqrt(1-e^2)*cos(Ea),0)
  P <- cos(w)*cos(Omega)-sin(w)*cos(i)*sin(Omega)
  Q <- sin(w)*cos(Omega)+cos(w)*cos(i)*sin(Omega)
  R <- cos(w)*sin(Omega)+sin(w)*cos(i)*cos(Omega)
  S <- cos(w)*cos(i)*cos(Omega)-sin(w)*sin(Omega)
  x <- o[,1]*P-o[,2]*Q
  y <- o[,1]*R+o[,2]*S
  z <- o[,1]*sin(w)*sin(i)+o[,2]*cos(w)*sin(i)
  vx <- vo[,1]*P-vo[,2]*Q
  vy <- vo[,1]*R+vo[,2]*S
  vz <- vo[,1]*sin(w)*sin(i)+vo[,2]*cos(w)*sin(i)
  cs <- cbind(x,y,z,vx,vy,vz)
  return(cs)
}

#Kepler to Cartesian
oe2cs <- function(oe,Mc=1){
  #deriv Ea by solving Kepler's equation M=E-e*sinE for E
  e <- oe[1]
  q <- oe[2]
  i <- oe[3]
  w <- oe[4]
  Omega <- oe[5]
  M <- oe[6]
  a <- q/(1-e)
  mu <- Mc*4*pi^2#6.67384e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2#m^3/s^2
  Ea <- kep.murison2(M,e)
  P <- c(cos(w)*cos(Omega)-sin(w)*cos(i)*sin(Omega),cos(w)*sin(Omega)+sin(w)*cos(i)*cos(Omega),sin(w)*sin(i))
  Q <- c(-sin(w)*cos(Omega)-cos(w)*cos(i)*sin(Omega),-sin(w)*sin(Omega)+cos(w)*cos(i)*cos(Omega),sin(i)*cos(w))
  R <- a*(cos(Ea)-e)*P + a*sqrt(1-e^2)*sin(Ea)*Q#position
  Edot <- sqrt(mu/a^3)/(1-e*cos(Ea))
  V <- -a*sin(Ea)*Edot*P + a*sqrt(1-e^2)*cos(Ea)*Edot*Q#velocity
  if(is.null(dim(R))){
      R <- matrix(R,nrow=1)
      V <- matrix(V,nrow=1)
  }
  cs <- cbind(R,V)
  return(cs)
}

###The following function have something wrong which make ele2cs and cs2oe not inversible
ele2cs <- function(oe){
  e <- oe[1]
  q <- oe[2]
  inc <- oe[3]
  w <- oe[4]
  Omega <- oe[5]
  M <- oe[6]
  a <- q/(1-e)
  mu <- 6.67384e-11*1.98855e30*1e-9*149597871^-3*(3600*24)^2*365.242199^2
  Ea <- kep.murison2(M,e)
  nu <- 2*atan(((1+e)/(1-e))^0.5*tan(Ea/2))
  p <- a*(1-e^2)
  r <- a*(1-e*cos(Ea))
  h <- (mu*a*(1-e^2))^0.5
  X <- r*(cos(Omega)*cos(w+nu)-sin(Omega)*sin(w+nu)*cos(inc))
  Y <- r*(sin(Omega)*cos(w+nu)+cos(Omega)*sin(w+nu)*cos(inc))
  Z <- r*sin(inc)*sin(w+nu)
  Vx <- X*h*e*sin(nu)/r/p-h/r*(cos(Omega)*sin(w+nu)+sin(Omega)*cos(w+nu)*cos(inc))
  Vy <- Y*h*e/r/p*sin(nu)-h/r*(sin(Omega)*sin(w+nu)-cos(Omega)*cos(w+nu)*cos(inc))
  Vz <- Z*h*e/r/p*sin(nu+h/r*sin(inc)*cos(w+nu))
  cs <- c(X,Y,Z,Vx,Vy,Vz)
  return(cs)
}
#funciton for deciding the new step size of LARKS model
h1 <- function(a1,r,Nstep){
  a0 <- 50000
  h0 <- a0^(3/2)/Nstep
  if(a1>0){
  h1 <- h0*(a1/a0)^(3/2)}else{
    h1 <- h0*(a0^(3/2)/(abs(a1)^(1/2)*r))}
}

Uutest <- function(u,U){
  val <- u[2]*U[1]-u[1]*U[2]-u[4]*U[3]+u[3]*U[4]
  if(val==0){
    print("The values of u and U are correct.")
  }else{
    print("The values of u and U are wrong!")}
}

orbit.HM <- function(element,Nstep,par.tide){
  Nper <- length(element[,1])-1
  time.step <- rep(NA,Nstep*Nper+1)
  phase <- array(data=NA,dim=c(Nper*Nstep+1,7),dimnames=c("Times","phase"))
  per.acc <- element[,1]
  e <- element[,2]
  q <- element[,3]
  I <- element[,4]
  w <- element[,5]
  Omega <- element[,6]
  a <- q/(1-e)
  period <- diff(per.acc)
  phase[1,1] <- 0
  phase[1,2:7] <- k2c(a[1],e[1],I[1],Omega[1],w[1],M=0,par.tide)
  time.step[1] <- 0
  for(i in 1:Nper){
    t <- seq(per.acc[i],per.acc[i+1],by=period[i]/Nstep)
    t <- t[-1]
    time.step[((i-1)*Nstep+2):(i*Nstep+1)] <- t
    n <- 2*pi/period[i]
    M <- n*(t-per.acc[i])
    for(j in 1:Nstep){
      phase[(i-1)*Nstep+j+1,2:7] <- k2c(a[i],e[i],I[i],Omega[i],w[i],M=M[j],par.tide)
    }
  }
  phase[,1] <- time.step
  E <- H0(phase[,2],phase[,3],phase[,4],phase[,5],phase[,6],phase[,7],par.tide)+H1(phase[,1],phase[,2],phase[,3],phase[,4],par.tide)
  return(list("phase"=phase,"E"=E))
}

oe2de <- function(oe,par.tide){
  de <- rep(NA,5)
  e <- oe[1]
  q <- oe[2]
  I <- oe[3]
  w <- oe[4]
  Omega <- oe[5]
  a <- q/(1-e)
  L <- sqrt(par.tide$mu*a)
  G <- L*sqrt(1-e^2)
  Theta <- G*cos(I)
  de <- c(w, Omega, L, G, Theta)
  return(de)
}

de2oe <- function(de, par.tide){
  w <- de[1]
  Omega <- de[2]
  L <- de[3]
  G <- de[4]
  Theta <- de[5]
  oe <- rep(NA,5)
  e <- sqrt(1-G^2/L^2)
  a <- L^2/par.tide$mu
  q <- a*(1-e)
  I <- acos(Theta/G)
  M <- 0
  oe <- c(e,q,I,w,Omega,M)
  return(oe)
}

ac <- function(par.comet,e){
  inda <- par.comet$inda
  inde <- par.comet$inde
  return(10^inda*(1-e)^inde)
}

#######This function is used to calcuate the orbit of minro bodies in the solar system perturbed by the Galactic tide.
orbit.integrator <- function(OE,tp,par.comet,par.tide,hmmethod,csmethod){
  Nstep <- par.comet$Nstep
  age <- par.comet$age
  qc <- par.comet$qc
  rc <- par.comet$rc
  ql <- par.comet$ql
  OE <- OE
  tp <- tp
  Nper <- par.comet$Nper
#  cat('Nper=',Nper,'\n')
  element <- c(tp,OE)
  Ntar <- 0
  Nesc <- 0
  Ttar <- NA
  Tesc <- NA
  as <- OE[2]/(1-OE[1])#element at last step
  qs <- OE[2]
  rs <- 2*as-qs
  de <- as.numeric(oe2de(OE,par.tide))
  DE <- c(w=de[1], Omega=de[2],L=de[3], G=de[4], Theta=de[5])
  for(k in 1:Nper){
      if(tp>age) break()
###OE ==> DE
      tper <- c(tp,tp+as^(3/2))#new times
      tper <- c(tp,tp + as^(3/2))
###integrator
      out <- as.numeric(ode(func=tide.HM,y=DE,parms=par.tide,times=tper,method=hmmethod)[2,])
      DE <- c(w=out[2], Omega=out[3],L=out[4], G=out[5], Theta=out[6])
###DE ==> OE
      OE <- de2oe(DE,par.tide)
      if(any(is.na(OE)) | any(is.na(DE))){
          cat('k=',k,'\n')
          cat('OE=',OE,'\n')
      }
      tp <- tper[2]
      element <- rbind(element,c(tp,OE))
  }
  return(list(element=element,Ttar=Ttar,Tesc=Tesc,Ntar=Ntar,Nesc=Nesc))
}

s2g <- function(x,y,z,b,l){
    xg <- cos(l)*cos(b)*x+sin(l)*y+cos(l)*sin(b)*z
    yg <- -sin(l)*cos(b)*x+cos(l)*y-sin(l)*sin(b)*z
    zg <- -sin(b)*x+cos(b)*z
    return(cbind(xg,yg,zg))
}

g2s <- function(x,y,z,b,l){
    xs <- cos(b)*cos(l)*x-cos(b)*sin(l)*y-sin(b)*z
    ys <- sin(l)*x+cos(l)*y
    zs <- cos(l)*sin(b)*x-sin(l)*sin(b)*y+cos(b)*z
    return(cbind(xs,ys,zs))
}

stellar2galactic <- function(dvs,l,b,a){
  bs <- -b
  if(l>=0 & l<=pi){
    ls <- l+pi
  }
  if(l>pi){
    ls <- l-pi
  }
  rotx <- matrix(data=c(1,0,0,0,cos(a),-sin(a),0,sin(a),cos(a)),nrow=3,ncol=3)
  roty <- matrix(data=c(cos(bs),0,sin(bs),0,1,0,-sin(bs),0,cos(bs)),nrow=3,ncol=3)
  rotz <- matrix(data=c(cos(ls),-sin(ls),0,sin(ls),cos(ls),0,0,0,1),nrow=3,ncol=3)
  dvg <- rotz%*%roty%*%rotx%*%dvs
  return(dvg)
}

galactic2stellar <- function(vg, ls, bs, a){
  rotx <- matrix(data=c(1,0,0,0,cos(a),sin(a),0,-sin(a),cos(a)),nrow=3,ncol=3)
  roty <- matrix(data=c(cos(bs),0,-sin(bs),0,1,0,sin(bs),0,cos(bs)),nrow=3,ncol=3)
  rotz <- matrix(data=c(cos(ls),sin(ls),0,-sin(ls),cos(ls),0,0,0,1),nrow=3,ncol=3)
  vs <- rotx%*%roty%*%rotz%*%vg
  return(vs)
}

####function to transfer orbital element in the ecliptic reference frame to Galactic coordinates
ele2bl <- function(oe){
   N <- as.integer(length(oe)/6)
   beta <- rep(NA,N)
   lambda <- rep(NA,N)
   for(i in 1:N){
       if(N==1){
           OE <- oe
       }else{
           OE <- oe[i,]
       }
      CS <- oe2cs(OE)
      bl <- xyz2bl.vec(CS[1],CS[2],CS[3])#its ecliptic coord.: lambda, beta
      betalambda <- ecliptic2gal(b=bl[1],l=bl[2])
#      bl <- equatorial2gal(alpha=betalambda[2],delta=betalambda[1])
#      bl <- equ2gal(alpha=betalambda[2],delta=betalambda[1])
      beta[i] <- betalambda[1]
      lambda[i] <- betalambda[2]
   }
   return(cbind(beta,lambda))
}

##transform orbital elements from ecliptic reference frame to the galactic reference frame
ele2gal <- function(oe){
   N <- as.integer(length(oe)/6)
   oe.gal <- oe
   for(i in 1:N){
       if(N==1){
           OE <- oe
       }else{
           OE <- oe[i,]
       }
       CS <- oe2cs(OE)
       Rnorm <- sqrt(CS[1:3]%*%CS[1:3])
       Vnorm <- sqrt(CS[4:6]%*%CS[4:6])
#       cat('CS=',CS,'\n')
       betalambda.pos <- xyz2bl.vec(CS[1],CS[2],CS[3])#its ecliptic coord.: lambda, beta
       betalambda.vel <- xyz2bl.vec(CS[4],CS[5],CS[6])#
       bl.pos <- ecliptic2gal(lambda=betalambda.pos[2],beta=betalambda.pos[1])
       bl.vel <- ecliptic2gal(lambda=betalambda.vel[2],beta=betalambda.vel[1])
       Rvec <- as.numeric(bl2xyz(bl.pos[1],bl.pos[2]))
       Vvec <-  as.numeric(bl2xyz(bl.vel[1],bl.vel[2]))
       CS.gal <- c(Rnorm*Rvec,Vnorm*Vvec)
       if(N==1){
           oe.gal <- cs2oe(CS.gal)
       }else{
           oe.gal[i,] <- cs2oe(CS.gal)
       }
   }
   return(oe.gal)
}

###transform orbital elements from galactic reference frame to ecliptic reference frame
ele2ecl <- function(oe){
   N <- as.integer(length(oe)/6)
   oe.ecl <- oe
   for(i in 1:N){
       if(N==1){
           OE <- oe
       }else{
           OE <- oe[i,]
       }
       CS <- oe2cs(OE)
       Rnorm <- sqrt(CS[1:3]%*%CS[1:3])
       Vnorm <- sqrt(CS[4:6]%*%CS[4:6])
       bl.pos <- xyz2bl.vec(CS[1],CS[2],CS[3])#its ecliptic coord.: lambda, beta
       bl.vel <- xyz2bl.vec(CS[4],CS[5],CS[6])#
       bl.pos <- gal2ecliptic(b=bl.pos[1],l=bl.pos[2])
       bl.vel <- gal2ecliptic(b=bl.vel[1],l=bl.vel[2])
       Rvec <- Rnorm*as.numeric(bl2xyz(bl.pos[1],bl.pos[2]))
       Vvec <- Vnorm*as.numeric(bl2xyz(bl.vel[1],bl.vel[2]))
       CS.ecl <- c(Rvec,Vvec)
       if(N==1){
           oe.ecl <- cs2oe(CS.ecl)
       }else{
           oe.ecl[i,] <- cs2oe(CS.ecl)
       }
   }
   return(oe.ecl)
}
####cartesian coordinates to longitude and latitude in spherical coordinates
xyz2bl <- function(x,y,z){
  b <- rep(NA,length(x))
  l <- rep(NA,length(x))
  for(j in 1:length(x)){
    l.test <- atan(y[j]/x[j])
    if(x[j]<0){
      l[j] <- l.test+pi
    }else if((x[j]>0) & (y[j]<0)){
      l[j] <- l.test+2*pi
    }else{
      l[j] <- l.test
    }
    b[j] <- atan(z[j]/sqrt(x[j]^2+y[j]^2))
  }
  return(cbind(b,l))
}

####Note that bl2xyz return coordinates in the heliocentric frame and the x axis point to the GC, yaxis point to the rotation of the Galaxy
bl2xyz <- function(b.rad,l.rad){
    x <- cos(b.rad)*cos(l.rad)
    y <- cos(b.rad)*sin(l.rad)
    z <- sin(b.rad)
    return(cbind(x,y,z))
}
###ref: http://physics.stackexchange.com/questions/88663/converting-between-galactic-and-ecliptic-coordinates
ecliptic2gal <- function(lambda,beta){
    rad2deg <- 180/pi
####epoch J2000
    lambdaG <- 180.01/rad2deg
    lambdaB <- 266.84/rad2deg
    betaG <- 29.80/rad2deg
    betaB <- -5.54/rad2deg
    BK <- 96.43/rad2deg
    z <- sin(betaG)*sin(beta)+cos(betaG)*cos(beta)*cos(lambda-lambdaG)
    y <- cos(beta)*sin(lambda-lambdaG)
    x <- cos(betaG)*sin(beta)-sin(betaG)*cos(beta)*cos(lambda-lambdaG)
    tmp <- xyz2bl.vec(x,y,z)
    bt <- tmp[1]
    lt <- (BK-tmp[2])%%(2*pi)
###test
    test <- FALSE
    if(test){
        cat('lambda=',lambda*rad2deg,'; beta=',beta*rad2deg,'\n')
        cat('l=',lt*rad2deg,'; b=',bt*rad2deg,'\n')
        zp <- sin(betaG)*sin(bt)+cos(betaG)*cos(bt)*cos(BK-lt)
        yp <- cos(bt)*sin(BK-lt)
        xp <- cos(betaG)*sin(bt)-sin(betaG)*cos(bt)*cos(BK-lt)
        tmp <- xyz2bl.vec(xp,yp,zp)
        betat <- tmp[1]
        lambdat <- tmp[2]+lambdaG
        cat('lambdat=',(lambdat*rad2deg)%%(360),'; betat=',betat*rad2deg,'\n')
    }
    return(c(bt,lt))
}


####galactic coordinates to ecliptic coordinates
gal2ecliptic <- function(b,l){
    rad2deg <- 180/pi
####epoch J2000
    lambdaG <- 180.01/rad2deg
    lambdaB <- 266.84/rad2deg
    betaG <- 29.80/rad2deg
    betaB <- -5.54/rad2deg
    BK <- 96.43/rad2deg
    zp <- sin(betaG)*sin(b)+cos(betaG)*cos(b)*cos(BK-l)
    yp <- cos(b)*sin(BK-l)
    xp <- cos(betaG)*sin(b)-sin(betaG)*cos(b)*cos(BK-l)
    tmp <- xyz2bl.vec(xp,yp,zp)
    beta <- tmp[,1]
    lambda <- tmp[,2]+lambdaG
    return(cbind(lambda,beta))
}

#####galactic coordinates to equatorial coordinates; the angle is in unit of rad.
gal2equ <- function(b,l){
####epoch 2000
    rad2deg <- 180/pi
    alphaG <- 192.85/rad2deg
    alphaB <- 266.40/rad2deg
    deltaG <- 27.13/rad2deg
    deltaB <- -28.94/rad2deg
    BK <- 122.9/rad2deg
    zp <- sin(deltaG)*sin(b)+cos(deltaG)*cos(b)*cos(BK-l)
    yp <- cos(b)*sin(BK-l)
    xp <- cos(deltaG)*sin(b)-sin(deltaG)*cos(b)*cos(BK-l)
    tmp <- xyz2bl.vec(xp,yp,zp)
    delta <- tmp[,1]
    alpha <- (tmp[,2]+alphaG)%%(2*pi)
    return(cbind(alpha,delta))
}

obs.lin.prop <- function(obs,t){
##t is in unit of day
##obs: ra(deg), dec(deg),plx(mas),pmra(mas/yr),pmdec(mas/yr)
    auyr2kms <- 4.74047
    kpcmyr2auyr <- 1e3*206265/1e6
    pc2au <- 206265
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    ra <- as.numeric(obs['ra'])/180*pi
    dec <- as.numeric(obs['dec'])/180*pi
    plx <- as.numeric(obs['parallax'])
    pmra <- as.numeric(obs['pmra'])
    pmdec <- as.numeric(obs['pmdec'])
    rv <- as.numeric(obs['radial_velocity'])

####propagation observables to states
    d <- as.numeric(1/plx)#kpc
    r <- as.numeric(bl2xyz(dec,ra))*d*1e3
    x <- r[1]
    y <- r[2]
    z <- r[3]
    vde <- pmdec*d
    vra <- pmra*d
    vp <- sqrt(vra^2+vde^2)
    vr <- rv/auyr2kms#au/yr
    vx.equ <- vr*cos(dec)*cos(ra)-vde*sin(dec)*cos(ra)-vra*sin(ra)##note: vr is positive if the star is moving away from the Sun
    vy.equ <- vr*cos(dec)*sin(ra)-vde*sin(dec)*sin(ra)+vra*cos(ra)
    vz.equ <- vr*sin(dec)+vde*cos(dec)
    x1 <- x+vx.equ*t/yr2d/pc2au
    y1 <- y+vy.equ*t/yr2d/pc2au
    z1 <- z+vz.equ*t/yr2d/pc2au

###convert time-varying states back to observables
    dec.ra <- xyz2bl.vec(x1,y1,z1)
    d1 <- sqrt(x1^2+y1^2+z1^2)*1e-3#kpc
    ra1.rad <- dec.ra[,2]#rad
    dec1.rad <- dec.ra[,1]
    ra1 <- ra1.rad*180/pi#deg
    dec1 <- dec1.rad*180/pi
###velocity to pm
    vequ <- array(NA,dim=c(length(t),3))
    for(j in 1:length(t)){
        rotz <- matrix(data=c(cos(ra1.rad[j]),sin(ra1.rad[j]),0.0,-sin(ra1.rad[j]),cos(ra1.rad[j]),0.0,0.0,0.0,1.0),nrow=3,ncol=3,byrow=TRUE)#o-xyz -> o-x'y'z'
                                        #roty <- matrix(data=c(cos(b.sun),0.0,sin(b.sun),0.0,1.0,0.0,-sin(b.sun),0.0,cos(b.sun)),nrow=3,ncol=3,byrow=TRUE)#o-x'y'z' ->o-x''y''z''
        roty <- matrix(data=c(cos(dec1.rad[j]),0.0,sin(dec1.rad[j]),0.0,1.0,0.0,-sin(dec1.rad[j]),0.0,cos(dec1.rad[j])),nrow=3,ncol=3,byrow=TRUE)
        vequ[j,] <- roty%*%rotz%*%as.numeric(c(vx.equ,vy.equ,vz.equ))
    }
    pmra1 <- vequ[,2]/d1#mas/yr
    pmdec1 <- vequ[,3]/d1#mas/yr
    rv1 <- vequ[,1]*auyr2kms
    out <- cbind(ra1,dec1,1/d1,pmra1,pmdec1,rv1)
    colnames(out) <- c('ra','dec','parallax','pmra','pmdec','radial_velocity')
    return(out)
}

#from observable to heliocentric xyzuvw
state.e2g <- function(state){
###from observables to xyzuvw with kpc and au/yr units
##input: gaia-like observables
##output: xyzuvw in units of kpc and au/yr
    auyr2kms <- 4.74047
    kpcmyr2auyr <- 1e3*206265/1e6
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    ra <- state[,'ra']/180*pi; dec <- state[,'dec']/180*pi; plx <- state[,'parallax']; pmra <- state[,'pmra']; pmdec <- state[,'pmdec']; rv <- state[,'radial_velocity']
    bl <- equ2gal(ra,dec)
    d <- 1/plx#kpc
    xyz <- bl2xyz(bl[,1],bl[,2])*d
###convert velocity
    vde <- pmdec*d/kpcmyr2auyr
    vra <- pmra*d/kpcmyr2auyr
    vp <- sqrt(vra^2+vde^2)
    vr <- rv/kpcmyr2kms#kpc/myr
    vx.equ <- vr*cos(dec)*cos(ra)-vde*sin(dec)*cos(ra)-vra*sin(ra)##note: vr is positive if the star is moving away from the Sun
    vy.equ <- vr*cos(dec)*sin(ra)-vde*sin(dec)*sin(ra)+vra*cos(ra)
    vz.equ <- vr*sin(dec)+vde*cos(dec)
    vhg <- e2g.vel(vx.equ,vy.equ,vz.equ)
    out <- cbind(xyz,vhg)
    colnames(out) <- c('x.kpc','y.kpc','z.kpc','u.kpcmyr','v.kpcmyr','w.kpcmyr')
    return(out)
}

#reliable function
state.g2e <- function(state,xyzuvw=FALSE){
    auyr2kms <- 4.74047
    kpcmyr2auyr <- 1e3*206265/1e6
    kpcmyr2kms <- kpcmyr2auyr*auyr2kms
    x <- state[,1]; y <- state[,2]; z <- state[,3]; vx <- state[,4]; vy <- state[,5]; vz <- state[,6]
    if(!xyzuvw){
        x <- -x
        y <- -y
        vx <- -vx
        vy <- -vy
    }
    N <- length(x)
    bl <- xyz2bl.vec(x,y,z)
    d <- sqrt(x^2+y^2+z^2)
    ad.rad <- gal2equ(bl[,1],bl[,2])
    ra.rad <- ad.rad[,1]
    dec.rad <- ad.rad[,2]
    ra <- ad.rad[,1]*180/pi
    dec <- ad.rad[,2]*180/pi
###velocity to pm
    bl <- xyz2bl.vec(vx,vy,vz)
    v <- sqrt(vx^2+vy^2+vz^2)
    vad <- gal2equ(bl[,1],bl[,2])
    vxe <- v*cos(vad[,2])*cos(vad[,1])
    vye <- v*cos(vad[,2])*sin(vad[,1])
    vze <- v*sin(vad[,2])
    vequ <- array(NA,dim=c(N,3))
    for(j in 1:N){
        rotz <- matrix(data=c(cos(ra.rad[j]),sin(ra.rad[j]),0.0,-sin(ra.rad[j]),cos(ra.rad[j]),0.0,0.0,0.0,1.0),nrow=3,ncol=3,byrow=TRUE)#o-xyz -> o-x'y'z'
                                        #roty <- matrix(data=c(cos(b.sun),0.0,sin(b.sun),0.0,1.0,0.0,-sin(b.sun),0.0,cos(b.sun)),nrow=3,ncol=3,byrow=TRUE)#o-x'y'z' ->o-x''y''z''
        roty <- matrix(data=c(cos(dec.rad[j]),0.0,sin(dec.rad[j]),0.0,1.0,0.0,-sin(dec.rad[j]),0.0,cos(dec.rad[j])),nrow=3,ncol=3,byrow=TRUE)
        vequ[j,] <- roty%*%rotz%*%c(vxe[j],vye[j],vze[j])
    }
    pm.ra <- vequ[,2]*kpcmyr2auyr/d#mas/yr
    pm.dec <- vequ[,3]*kpcmyr2auyr/d#mas/yr
    rv <- vequ[,1]*kpcmyr2kms#km/s
    out <- cbind(ra,dec,1/d,pm.ra,pm.dec,rv)
    colnames(out) <- c('ra','dec','parallax','pmra','pmdec','radial_velocity')
    return(out)
}

#####Equatorial coordinates to galactic coordinates; the angle is in unit of rad.
equ2gal <- function(alpha,delta){
    rad2deg <- 180/pi
    alphaG <- 192.85/rad2deg
    alphaB <- 266.40/rad2deg
    deltaG <- 27.13/rad2deg
    deltaB <- -28.94/rad2deg
    BK <- 122.9/rad2deg
    z <- sin(deltaG)*sin(delta)+cos(deltaG)*cos(delta)*cos(alpha-alphaG)
    y <- cos(delta)*sin(alpha-alphaG)
    x <- cos(deltaG)*sin(delta)-sin(deltaG)*cos(delta)*cos(alpha-alphaG)
    tmp <- xyz2bl.vec(x,y,z)
    bt <- tmp[,1]
    lt <- (BK-tmp[,2])%%(2*pi)
    test <- FALSE
###test
    if(test){
        cat('alpha=',alpha*rad2deg,'; delta=',delta*rad2deg,'\n')
        cat('b=',bt*rad2deg,'; l=',lt*rad2deg,'\n')
        zp <- sin(deltaG)*sin(bt)+cos(deltaG)*cos(bt)*cos(BK-lt)
        yp <- cos(bt)*sin(BK-lt)
        xp <- cos(deltaG)*sin(bt)-sin(deltaG)*cos(bt)*cos(BK-lt)
        tmp <- xyz2bl.vec(xp,yp,zp)
        deltat <- tmp[1]
        alphat <- (tmp[2]+alphaG)%%(2*pi)
        cat('alphat=',(alphat*rad2deg)%%360,'; deltat=',deltat*rad2deg,'\n')
    }
    return(cbind(bt,lt))
}



###ref: http://www2.astro.psu.edu/users/rbc/a501/c1_spherical_astronomy.pdf
equtorial2gal <- function(alpha,delta){
####epoch J1950
    alpha0 <- 282.25/rad2deg
    delta0 <- 62.6/rad2deg
    l0 <- 33/rad2deg
    x <- cos(delta)*cos(alpha-alpha0)
    y <- cos(delta)*sin(alpha-alpha0)*cos(delta0)+sin(delta)*sin(delta0)
    z <- sin(delta)*cos(delta0)-cos(delta)*sin(alpha-alpha0)*sin(delta0)
    tmp <- xyz2bl(x,y,z)
    bt <- tmp[1]
    lt <- (tmp[2]+l0)%%(2*pi)
    return(c(bt,lt))
}

####equatorial velocity to velocity in galactic coordinates in the helio-static frame
e2g.vel <- function(equ){
    if(is.null(dim(equ))) equ <- t(equ)
    x <- equ[,1]
    y <- equ[,2]
    z <- equ[,3]
    N <- length(x)
    r <- sqrt(x^2+y^2+z^2)
    ad <- xyz2equ(x,y,z)
    bl <- equ2gal(ad[,1],ad[,2])
    rgal <- r*bl2xyz(bl[,1],bl[,2])
    return(rgal)
}

####vector in galactic velocity to equatorial velocity in the helio-centric frame
g2e.vel <- function(gal){
    if(is.null(dim(gal))) gal <- t(gal)
    x <- gal[,1]
    y <- gal[,2]
    z <- gal[,3]
    r <- sqrt(x^2+y^2+z^2)
    bl <- xyz2bl.vec(x,y,z)
    ad <- gal2equ(bl[,1],bl[,2])
    vxe <- r*cos(ad[,2])*cos(ad[,1])
    vye <- r*cos(ad[,2])*sin(ad[,1])
    vze <- r*sin(ad[,2])
    return(cbind(vxe,vye,vze))
}

#####velocity kick from stellar encounters using classical impulse approximation (Rickman 1976)
dv.gal <- function(b.target,b.sun,mstar,vph){
###only allow one target, but could be multiple encounters
###b.target, b.sun are in unit of au, mstar is in unit of Msun, vph is in unit of au/yr
    N <- round(length(b.target)/3)
    if(N==1){
        b.target <- matrix(b.target,nrow=1)
        b.sun <- matrix(b.sun,nrow=1)
    }
    bt2 <- b.target[,1]^2+b.target[,2]^2+b.target[,3]^2
    bs2 <- b.sun[,1]^2+b.sun[,2]^2+b.sun[,3]^2
    dvx <- 2*mu/vph*(b.target[,1]/bt2-b.sun[,1]/bs2)
    dvy <- 2*mu/vph*(b.target[,2]/bt2-b.sun[,2]/bs2)
    dvz <- 2*mu/vph*(b.target[,3]/bt2-b.sun[,3]/bs2)
#    cat('length(dvx)=',length(dvx),'\n')
    return(cbind(dvx,dvy,dvz))
}
deg2hdms <- function(RAdeg,DEdeg){
    val <- RAdeg/360*24#hour
    RAh <- floor(val)
    RAm <- floor((val%%1)*60)
    RAs <- (((val%%1)*60)%%1)*60
    sig <- sign(DEdeg)
    DEdeg <- abs(DEdeg)
    DEd <- sig*floor(DEdeg)
    DEm <- floor((DEdeg%%1)*60)
    DEs <- (((DEdeg%%1)*60)%%1)*60
    return(cbind(RAh,RAm,RAs,DEd,DEm,DEs))
}

hdms2deg <- function(rah,ram,ras,ded,dem,des,sn=NULL){
    RA.deg <- (rah+ram/60+ras/3600)/24*360
    if(is.null(sn)){
        sn <- sign(ded)
    }
    DEC.deg <- sn*(abs(ded)+dem/60+des/3600)
    RA.rad <- RA.deg/180*pi
    DEC.rad <- DEC.deg/180*pi
    cbind(RA.deg,DEC.deg,RA.rad,DEC.rad)
}

###transform from mu, PA to mu.ra and mu.dec
mu2pm <- function(mu.masyr, PA.deg){
    mu.ra <- mu.masyr*cos(PA.deg/180*pi)
    mu.dec <- mu.masyr*sin(PA.deg/180*pi)
    cbind(mu.ra,mu.dec)
}
####error transformation
emu2pm <- function(mu,PA,emu, ePA){
    theta <- PA/180*pi
    etheta <- ePA/180*pi
    emu.ra <- sqrt((emu*cos(theta))^2+(mu*sin(theta)*etheta)^2)
    emu.dec <- sqrt((emu*sin(theta))^2+(mu*cos(theta)*etheta)^2)
    cbind(emu.ra,emu.dec)
}

lp <- function(plx,pmra,pmdec,rv){
    pc2au <- 206265
    c <- 4.74047
    mu <- sqrt(pmra^2+pmdec^2)
    vt <- c*mu/plx#km/s
    v <- sqrt(vt^2+rv^2)
    sin.theta <- vt/v
    cos.theta <- -rv/v
    dph <- 1000/plx*sin.theta#pc
    tph <- pc2au*1000/plx*cos.theta/(v/c)/1e6#Myr
    cbind(dph,tph)
}

lp.pair <- function(X1,Y1,Z1,vx1,vy1,vz1,X2,Y2,Z2,vx2,vy2,vz2){
    x <- X2-X1
    y <- Y2-Y1
    z <- Z2-Z1
    vx <- vx2-vx1
    vy <- vy2-vy1
    vz <- vz2-vz1
    tph <- -(x*vx+y*vy+z*vz)/(vx^2+vy^2+vz^2)
    xp <- x+vx*tph
    yp <- y+vy*tph
    zp <- z+vz*tph
    dph <- sqrt(xp^2+yp^2+zp^2)
    cbind(dph,tph,xp,yp,zp,vx,vy,vz)
}

lp <- function(t,r,v){
    if(is.null(dim(r))){
        r <- matrix(r,nrow=1)
        v <- matrix(v,nrow=1)
    }
    tph <- -rowSums(r*v)/rowSums(v^2)
    rph <- r+v*tph
    dph <- sqrt(rowSums(rph^2))
    cbind(dph,t+tph,rph,v)
}

rot2gal <- function(R,Rdot,phi,phidot,z,zdot){
    x <- R*cos(phi)
    y <- R*sin(phi)
    vx <- Rdot*cos(phi)-R*phidot*sin(phi)
    vy <- Rdot*sin(phi)+R*phidot*cos(phi)
    cbind(x,y,z,vx,vy,zdot)
}
tcol <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1,], rgb.val[2,], rgb.val[3,],
                 max = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    invisible(t.col)
}
calc.kick <- function(xs,ys,zs,vxs,vys,vzs,Ms,ind.bary,enc,verbose=FALSE){
    xb <- xs[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
    yb <- ys[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
    zb <- zs[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
    vxb <- vxs[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
    vyb <- vys[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
    vzb <- vzs[ind.bary]%*%Ms[ind.bary]/sum(Ms[ind.bary])
    xsb <- as.numeric(xs-xb)
    ysb <- as.numeric(ys-yb)
    zsb <- as.numeric(zs-zb)
    vxsb <- as.numeric(vxs-vxb)
    vysb <- as.numeric(vys-vyb)
    vzsb <- as.numeric(vzs-vzb)
    xes <- outer(enc[,2],xsb,'-')
    yes <- outer(enc[,3],ysb,'-')
    zes <- outer(enc[,4],zsb,'-')
    vxes <- outer(enc[,5],xsb,'-')
    vyes <- outer(enc[,6],ysb,'-')
    vzes <- outer(enc[,7],zsb,'-')
    Des <- cbind(xes,yes,zes)
    Ves <- cbind(vxes,vyes,vzes)
    des <- sqrt(xes^2+yes^2+zes^2)
    ves <- sqrt(vxes^2+vyes^2+vzes^2)
    renc <- enc[,2:4]
    mu <- (2*pi)^2*enc[,'Menc']
    if(LA & sim.type=='tide' & ind11==0){
        Psi <- acos(1/sqrt(1+(des*kpc2au)^2*(ves*kpcmyr2auyr)^4/mu^2))
        phi <- 2*Psi-pi
        dvx <- sum(vxes*cos(phi)-ves*xes/des*sin(phi)-vxes)
        dvy <- sum(vyes*cos(phi)-ves*yes/des*sin(phi)-vyes)
        dvz <- sum(vzes*cos(phi)-ves*zes/des*sin(phi)-vzes)
    }else{
        a <- 2*mu/(ves*kpcmyr2auyr)#2*G*M/venc
                                        #    dvx <- a%*%(xes/des^2-enc[,2]/enc[,11]^2)/kpcmyr2auyr/kpc2au#kpc/myr
        dvx <- colSums(a*(xes/des^2)/kpcmyr2auyr/kpc2au)
                                        #    dvy <- a%*%(yes/des^2-enc[,3]/enc[,11]^2)/kpcmyr2auyr/kpc2au#kpc/myr
        dvy <- colSums(a*(yes/des^2)/kpcmyr2auyr/kpc2au)
                                        #    dvz <- a%*%(zes/des^2-enc[,4]/enc[,11]^2)/kpcmyr2auyr/kpc2au#kpc/myr
        dvz <- colSums(a*(zes/des^2)/kpcmyr2auyr/kpc2au)
    }
    tmp <- t(rbind(dvx,dvy,dvz))
#    cat('dv=',sqrt(dvx^2+dvy^2+dvz^2)*kpcmyr2kms,'km/s\n')
    if(verbose){
        cat('dvx=',dvx*kpcmyr2kms,'km/s\n')
        cat('dvy=',dvy*kpcmyr2kms,'km/s\n')
        cat('dvz=',dvz*kpcmyr2kms,'km/s\n')
        cat('denc=',enc[,11]*1e3,'pc\n')
        cat('des=',range(des)*1e3,'pc\n')
        cat('dsb=',sqrt(xsb^2+ysb^2+zsb^2)*1e3,'pc\n')
        cat('venc=',enc[,10]*kpcmyr2kms,'\n')
        cat('range(tmp)=',range(tmp)*kpcmyr2kms,'km/s\n\n')
    }
    return(tmp)
}
ra.dec.conv <- function(ra,dec,pmra,pmdec,epoch.from,epoch.to=2000){
#convert to epoch 2000; ra, dec in deg; pmra, pmdec in mas/yr
    depoch <- epoch.from-epoch.to
    dra <- pmra/cos(dec/180*pi)/3600000*depoch
    ddec <- pmdec/3600000*5
    ra <- ra-dra
    dec <- dec-ddec
    return(cbind(ra,dec))
}

u2cov <- function(ut,err,u){
    cors <- covs <- array(NA,dim=c(nrow(ut),10))
    for(j in 1:nrow(ut)){
        c1 <- as.numeric(c(ut[j,1],rep(0,4)))
        c2 <- as.numeric(c(ut[j,2:3],rep(0,3)))
        c3 <- as.numeric(c(ut[j,4:6],rep(0,2)))
        c4 <- as.numeric(c(ut[j,7:10],rep(0,1)))
        c5 <- as.numeric(ut[j,11:15])
        U <- cbind(c1,c2,c3,c4,c5)
        cov.mat <- solve(t(U)%*%U)
        if(u[j]>1) cov.mat <- cov.mat*u[j]^2
        cor.mat <- cov.mat/outer(as.numeric(err[j,]),as.numeric(err[j,]),"*")
        covs[j,] <- as.numeric(c(cov.mat[1,2:5],cov.mat[2,3:5],cov.mat[3,4:5],cov.mat[4,5]))
        cors[j,] <- as.numeric(c(cor.mat[1,2:5],cor.mat[2,3:5],cor.mat[3,4:5],cor.mat[4,5]))
    }
    return(list(cov=covs,cor=cors))
}

delete.overlap <- function(ra,dec,cat){
    ind.keep <- c()
    ind.rm <- c()
    for(j in 1:nrow(tab)){
        if(j%%1000==0) cat('j=',j,'\n')
        if(!any(j==ind.rm)){
            inds <- (1:nrow(tab))[-j]
            index <- which(abs(ra[-j]-ra[j])<5/3600 & abs(dec[-j]-dec[j])<5/3600)
            if(length(index)>0){
                cat('j=',j,'\n')
                cat('length(index)=',length(index),'\n')
                ii <- c(j,inds[index])
                tmp <- cat[ii]
                pp <- which.min(match(tmp,cats))
                ind.keep <- c(ind.keep,ii[pp])
                ind.rm <- c(ind.rm,ii[-pp])
                cat('cat:',cat[ii[pp]],'\n')
            }
        }
    }
    return(tab)
}

bs.step <- function(f,time,y0,Dt,parameters,tol=1e-7){
##t are forward times
    if(all(time<=0)){
        sign <- -1
    }
    if(all(time>=0)){
        sign <- 1
    }
    tt <- abs(time)
    tmp <- array(NA,dim=c(length(tt),length(y0)))
    Nstep <- ceiling((max(tt)-min(tt))/Dt)
#    if(mean(diff(time))==Dt) Nstep <- Nstep-1
    ic <- y0
    for(j in 1:Nstep){
        inds <- which(tt>=(j-1)*Dt & tt<=(j*Dt))
        pp <- bulirsch_stoer(f=f, t=sign*(tt[inds]-min(tt[inds])), y0=ic,tol=tol,parameters=parameters)
        ic <- t(pp[nrow(pp),,drop=FALSE])
        tmp[inds,] <- pp
    }
    return(tmp)
}

binG <- function(ts,gs,dN,method='sum'){
    inds <- sort(ts,index.return=TRUE)$ix
    ts <- ts[inds]
    gs <- gs[inds]
    Nbin <- ceiling(length(gs)/dN)
    gg <- rep(NA,Nbin)
    tt <- rep(NA,Nbin)
    for(k in 1:Nbin){
        inds <- ((k-1)*dN+1):min(k*dN,length(ts))
        if(method=='sum'){
            gg[k] <- sum(gs[inds])
        }else if(method=='max'){
            gg[k] <- max(gs[inds])
        }
        tt[k] <- median(ts[inds])
    }
    cbind(tt,gg)
}

error.ellipse <- function(x,y,covar,percent=95){
#http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
#http://www.r-tutor.com/elementary-statistics/probability-distributions/chi-squared-distribution
    tmp <- eigen(covar)
    lambda <- abs(tmp$values)
    v1 <- tmp$vectors[,1]
    angle <- atan(v1[2]/v1[1])
    if(angle < 0) angle <- angle+2*pi
    s <- qchisq(percent/100, df=2)
    a <- sqrt(s*lambda[1])
    b <- sqrt(s*lambda[2])
    t <- seq(0, 2*pi, by=pi/1000)
    xc <- x
    yc <- y
    xt <- a*cos(t)
    yt <- b*sin(t)
    rot <- array(c(cos(angle),sin(angle),-sin(angle),cos(angle)),dim=c(2,2))
    xy <- rot%*%rbind(xt,yt)
    xt <- xy[1,]+xc
    yt <- xy[2,]+yc
    return(cbind(xt,yt))
}
ellipse.edge <- function(covar,xs,ys,ic){
                                        # covariance of x,y
                                        # xs and ys are the coordinates of the center of the ellipse and two nearby points on the reference orbit to determine the slope
                                        # ... are any arguments that can be passed to function lines
    tmp <- eigen(covar)
    lambda <- abs(tmp$values)
    v1 <- tmp$vectors[,1]
    angle <- atan(v1[2]/v1[1])
    if(angle < 0) angle <- angle+2*pi
    a <- sqrt(5.991*lambda[1])
    b <- sqrt(5.991*lambda[2])
    t <- seq(0, 2*pi, by=pi/1000)
    xc <- xs[ic]
    yc <- ys[ic]
    xt <- a*cos(t)
    yt <- b*sin(t)
    rot <- array(c(cos(angle),sin(angle),-sin(angle),cos(angle)),dim=c(2,2))
    xy <- rot%*%rbind(xt,yt)
    xt <- xy[1,]+xc
    yt <- xy[2,]+yc
    s1 <- mean(diff(ys))/mean(diff(xs))
    ind.up <- which(yt>yc)
    ind.low <- which(yt<=yc)
    s2up <- diff(yt[ind.up])/diff(xt[ind.up])
    s2low <- diff(yt[ind.low])/diff(xt[ind.low])
    ind1 <- ind.up[which.min(abs(s2up-s1))]
    ind2 <- ind.low[which.min(abs(s2low-s1))]
    xup <- xt[ind1]
    xlow <- xt[ind2]
    yup <- yt[ind1]
    ylow <- yt[ind2]
###determine semi-ellipse divided by the orthogonal line
    slope <- (yup-ylow)/(xup-xlow)
    intercept <- (ylow*xup-xlow*yup)/(xup-xlow)
    ind.semi.up <- which(yt>=(intercept+slope*xt))
    ind.semi.low <- which(yt<=(intercept+slope*xt))
    xy <- cbind(xt[c(ind1,ind2)],yt[c(ind1,ind2)])
###sort the index
    t2 <- (t+angle)%%(2*pi)
    ind.semi.up <- ind.semi.up[sort(t2[ind.semi.up],index.return=TRUE)$ix]
    ind.semi.low <- ind.semi.low[sort(t2[ind.semi.low],index.return=TRUE)$ix]
    xs.up <- xt[ind.semi.up]
    xs.low <- xt[ind.semi.low]
    ys.up <- yt[ind.semi.up]
    ys.low <- yt[ind.semi.low]
##output
    return(list(edge=xy,xt=xt,yt=yt,xup=xs.up,yup=ys.up,xlow=xs.low,ylow=ys.low))
#    lines(xt,yt,...)
#    points(xt[c(ind1,ind2)],yt[c(ind1,ind2)],...)
}
rm.out <- function(x,y,alpha=3){
    r <- sqrt(x^2+y^2)
    ind.out <- which(abs(diff(r))>abs(median(diff(r)))+alpha*sd(diff(r)))
    ind.in <- (1:length(r))[-ind.out]
    for(k in ind.out){
        ind.rep <- ind.in[which.min(abs(ind.in-k))]
        x[k] <- x[ind.rep]
        y[k] <- y[ind.rep]
    }
    return(cbind(x,y))
}

#calPCA function(cov.mat){
#    Ntime <- dim(cov.mat)[3]
#    sapply(1:Ntime,function(i) sqrt(sum(eigen(cov.out[1:3, 1:3, i])$values)))
#}
