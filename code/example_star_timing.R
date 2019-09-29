###test
if(Par$star=='XO3'){
    pdf('XO3.pdf',6,6)
    par(mar=c(5,5,1,1))
                                        #par(mfrow=c(2,2))
    tab <- read.table('../input/XO3_transit_W08.dat')
    t <- tab[,1]
    ddt <- (tab[,1]-tab[2,1])%%KeplerPar['pb']
    dt[dt>KeplerPar['pb']/2] <- dt[dt>KeplerPar['pb']/2]-KeplerPar['pb']
    edt <- tab[,2]
    y <-OutTime$out$Tc[,2]-OutTime$out$Tc[1,2]
    ##Tc
                                        #plot(jd.utc,y,xlab='jd[utc]',ylab='dTc')
                                        #plot(jd.utc[-1],diff(OutTime$out$Tc[,2]),xlab='jd[utc]',ylab='dTc')
    y <- y*DAYSEC
    ddt <- ddt*DAYSEC
    edt <- edt*DAYSEC
    t2 <- rowSums(bjd.tdb)
    plot(t2,y,xlab=expression(BJD[TDB]),ylab=expression(Delta*T[c]*' [s]'),ylim=range(dt,y),type='l')
    lines(t2,y*1e2,xlab=expression(BJD[TDB]),lty=2)
    points(t,dt,col='red')
    arrows(t,dt-edt,t,dt+edt,length=0.05,angle=90,code=3,col='red')
    text(x=2454650,y=20,labels=expression(Delta*T[c]^{GR}))
    text(x=2454650,y=-60,labels=expression(100*Delta*T[c]^{GR}))
    dev.off()
}

if(Par$star=='XO3'){
    vp <- vTC-rowSums(vTC*uOT)*uOT
    v <- sqrt(rowSums(vp^2))*auyr2kms
#    tab <- read.table('../input/Kepler210c_TDV_TPV.dat',header=TRUE)
    pdf('XO3_TDV.pdf',6.5,4)
    par(mar=c(5,5,1.1,1.1),oma=c(0,0,0,12))
#    R <- as.numeric(par.comp['Rstar'])*(1+2*as.numeric(par.comp['Rratio']))*Rsun.km
    dur <- Rsun.km/v/3600#hour
    dur.frac <- (dur-mean(dur))/mean(dur)
    tt <- rowSums(bjd.tdb)
    dur2 <- Rsun.km/(OutTime$out$Vc*auyr2kms)/3600
    dur2.frac <- (dur2-mean(dur2))/mean(dur2)
    ylab <- 'Transit Duration Fraction'
    t0 <- 2455000
    plot(tt-t0,dur.frac,xlab=expression(BJD[TDB]-2455000),ylab=ylab,ylim=ylim,type='l',col='red')
###superpose model prediction
    lines(tt-t0,dur2.frac,col='blue')
    lines(tt-t0,dur2.frac*1e2,col='red')
    abline(h=0,lty=2)
#    legend('topright',inset=c(-0.9,0),legend=c('Raw TDF data','Binned TDF data','Relativistic TDF','Amplified relativistic TDF'),col=c('grey','black','blue','red'),bty='n',xpd=NA,pch=c(1,20,NA,NA),lty=1)
    dev.off()
}

if(Par$star=='kepler210'){
#    v <- sqrt(rowSums(vTC^2))*auyr2kms
    vp <- vTC-rowSums(vTC*uOT)*uOT
    v <- sqrt(rowSums(vp^2))*auyr2kms
    tab <- read.table('../input/Kepler210c_TDV_TPV.dat',header=TRUE)
#    pdf('kepler210_TDV.pdf',9,4)
#    par(mar=c(5,5,1.1,1.1),mfrow=c(1,2),oma=c(0,0,0,8))
    pdf('kepler210_TDV.pdf',6.5,4)
    par(mar=c(5,5,1.1,1.1),oma=c(0,0,0,12))
    R <- as.numeric(par.comp['Rstar'])*(1+2*as.numeric(par.comp['Rratio']))*Rsun.km
    dur <- R/v/3600#hour
    dur.frac <- (dur-mean(dur))/mean(dur)
    tt <- rowSums(bjd.tdb)
    dur2 <- R/(OutTime$out$Vc*auyr2kms)/3600
    dur2.frac <- (dur2-mean(dur2))/mean(dur2)
    for(j in 1){
        if(j==1){
            ylim <- range(dur.frac,tab[,2])
        }else{
            ylim <- c(-0.002,0.002)
        }
        ylab <- 'Transit Duration Fraction'
        t0 <- 2455000
    plot(tt-t0,dur.frac,xlab=expression(BJD[TDB]-2455000),ylab=ylab,ylim=ylim,type='l',col='red')
    points(tab[,1]-t0,tab[,2],col='grey')
    arrows(tab[,1]-t0,tab[,2]-tab[,3],tab[,1]-t0,tab[,2]+tab[,3],length=0.05,angle=90,code=3,col='grey')
###binning
    data <- wtb(tab[,1],tab[,2],tab[,3],dt=100,sj=0)
    points(data[,1]-t0,data[,2],col='black',pch=20)
    arrows(data[,1]-t0,data[,2]-data[,3],data[,1]-t0,data[,2]+data[,3],length=0.05,angle=90,code=3,col='black',pch=20)
###superpose model prediction
        lines(tt-t0,dur2.frac,col='blue')
        lines(tt-t0,dur2.frac*1e2,col='red')
#    lines(tt-t0,dur.frac,col='red')
#    abline(h=0,lty=2)
}
    legend('topright',inset=c(-0.9,0),legend=c('Raw TDF data','Binned TDF data','Relativistic TDF','Amplified relativistic TDF'),col=c('grey','black','blue','red'),bty='n',xpd=NA,pch=c(1,20,NA,NA),lty=1)
    dev.off()
}
