dir.out <- '../results/'
if(!file.exists(dir.out)) system(paste('mkdir',dir.out))
#plot.type <- 'all'
plot.type <- 'useful'
jd.utc <- rowSums(utc)
Ntime <- nrow(utc)
#source('compare_with_T2.R')

dt <- diff(rowSums(utc))[1]
tS <- OutTime$tS
fname <- paste0(Par$star,'_timing_',Par$BinaryModel,'_Ntime',Par$Nepoch,'.pdf')
if(Par$binary){
    fname <- paste0(Par$star,'_timing_',Par$BinaryModel,'_dt',round(dt),'day_Ntime',Par$Nepoch,'_Omega',round(Par$Omega*180/pi),'_omegaT',round(Par$omegaT*180/pi),'.pdf')
}
Dt <-rowSums(OutTime$tauE-tS)#day
utcs <- rowSums(utc)
tcbs <-rowSums(OutTime$BJDtcb)
#t1 <- (tdbs-min(tdbs))/DJY
t1 <- tcbs-Par$tpos

fout <- paste0(dir.out,'AllTimes_',fname)
cat(fout,'\n')
pdf(fout,16,16)
#
size <- 1.5
par(cex=size,cex.lab=size,cex.axis=size,mar=c(6,6,3,1),cex.main=size)
layout(matrix(data=c(rep(1,4),rep(2,4),rep(3,4),rep(4,4)),nrow=4,ncol=4,byrow = TRUE))

#t <- jd.utc-min(jd.utc)
t <- time_Jd2yr(utc)
#xlab=expression(JD[UTC]-JD[UTC0])
ylabs <- c('UTC-TAI','TCG-TAI','TT-TAI','UT1-TAI','TDB-TAI','TCB-TAI','BJD[TDB]-TAI','BJD[TCB]-TAI','tB-TAI','tauE-TAI')
tai <- OutObs$JDtai
dutc <- time_T2mT2(utc,tai)
dtcg <- time_T2mT2(OutObs$JDtcg,tai)
dtt <- time_T2mT2(OutObs$JDtt,tai)
dut1 <- time_T2mT2(OutObs$JDut1,tai)
dtdb <- time_T2mT2(OutObs$JDtdb,tai)
dtcb <- time_T2mT2(OutObs$JDtcb,tai)
dBJDtdb <- time_T2mT2(OutTime$BJDtdb,tai)
dBJDtcb <- time_T2mT2(OutTime$BJDtcb,tai)
dtB <- time_T2mT2(OutTime$tB,tai)
dtauE <- time_T2mT2(OutTime$tauE,tai)
xlim <- range(min(t),max(t)+0.2*(max(t)-min(t)))
xlab <- 'JY[UTC] (year)'
cols <- c('black','blue','red','orange','green','brown','cyan','steelblue','pink','darkgrey','yellow')
for(j in 1:4){
    if(j==1){
        ylim <- range(dtt,dtdb)
        tit <- 'Zoom-in for TT and TDB'
    }else if(j==2){
        ylim <- range(dutc,dtcg,dtt,dut1,dtdb,dtcb)
        tit <- 'Zoom-in for UTC, TCG, TT, UT1, TDB and TCB'
    }else if(j==3){
        ylim <- range(dutc,dtcg,dtt,dut1,dtdb,dtcb,dBJDtcb,dBJDtdb)
        tit <- 'Zoom-in for UTC, TCG, TT, UT1, TDB, TCB, BJDtcb and BJDtdb'
    }else{
        ylim <- range(dutc,dtcg,dtt,dut1,dtdb,dtcb,dBJDtcb,dBJDtdb,dtB,dtauE)
        tit <- 'For all times'
    }
plot(t,dutc,xlab=xlab,ylab='Time - TAI [s]',type='l',ylim=ylim,xlim=xlim,col=cols[1],main=tit)
abline(h=0,col='black',lty=3,lwd=3)
abline(v=1960,col='grey',lty=2)
lines(t,dtcg,col=cols[2])
lines(t,dtt,col=cols[3])
lines(t,dut1,col=cols[4])
lines(t,dtdb,col=cols[5])
lines(t,dtcb,col=cols[6])
lines(t,dBJDtdb,col=cols[7])
lines(t,dBJDtcb,col=cols[8])
lines(t,dtB,col=cols[9])
lines(t,dtauE,col=cols[10])
legend('topright',legend=ylabs,col=cols,lty=1,bty='n',cex=size)
}

par(mar=c(5,5,1,1),mfrow=c(4,4))
plot(t,OutTime$TropoDelay,xlab=xlab,ylab='Tropospheric delay [s]',type='l')
plot(t,OutTime$RoemerSolar,xlab=xlab,ylab='RoemerSolar [s]',type='l')
plot(t,time_T2mT2(OutObs$JDtcb,OutObs$JDtt),xlab=xlab,ylab='EinsteinSolar [s]',type='l')
plot(t,OutTime$ShapiroSolar,xlab=xlab,ylab='ShapiroSolar [s]',type='l')
plot(t,OutTime$ShapiroPlanet$Sun,xlab=xlab,ylab='ShapiroSun [s]',type='l')
plot(t,OutTime$ShapiroPlanet$Earth,xlab=xlab,ylab='ShapiroEarth [s]',type='l')
plot(t,OutTime$ShapiroPlanet$Saturn,xlab=xlab,ylab='ShapiroSaturn [s]',type='l')
plot(t,OutTime$ShapiroPlanet$Jupiter,xlab=xlab,ylab='ShapiroJupiter [s]',type='l')
if(Par$binary & Par$Np>0){
    if(any(names(Par)=='T0')){
        dt <- (rowSums(OutTime$tB)-Par$T0)/DJY
    }
    if(any(names(Par)=='Tp')){
        dt <- (rowSums(OutTime$tB)-Par$Tp)/DJY
    }
    if(any(names(Par)=='Tasc')){
        if(!is.na(Par$Tasc)){
            dt <- (rowSums(OutTime$tB)-Par$Tasc)/DJY
        }
    }
    P <- exp(Par$logP)
    phase <- (dt%%P)/P
    ind <-sort(phase,index.return=TRUE)$ix
    plot(t,OutTime$RoemerTarget,xlab=xlab,ylab='RoemerTarget [s]',type='l')
    plot(phase[ind],OutTime$RoemerTarget[ind],xlab='Phase',ylab='RoemerTarget [s]',type='l')
    plot(t,OutTime$EinsteinTarget,xlab=xlab,ylab='EinsteinTarget [s]',type='l')
    plot(phase[ind],OutTime$EinsteinTarget[ind],xlab='Phase',ylab='EinsteinTarget [s]',type='l')
    plot(t,OutTime$ShapiroTarget,xlab=xlab,ylab='ShapiroTarget [s]',type='l')
#    use tB coordinate time as the x axis
    if(any(names(Par)=='Tp')){
        t0 <- as.numeric(time_Jd2yr(cbind(Par$Tp,0)))
    }else if(any(names(Par)=='T0')){
        t0 <- as.numeric(time_Jd2yr(cbind(Par$T0,0)))
    }else if(any(names(Par)=='Tasc')){
        t0 <- as.numeric(time_Jd2yr(cbind(Par$Tasc,0)))
    }
    plot(phase[ind],OutTime$ShapiroTarget[ind],xlab='Phase',ylab='ShapiroTarget [s]',type='l')
    if(Par$star=='PSRJ0740+6620' & FALSE){
        tmp0 <- readLines('../input/J0740+6620_NANOGrav_12yv2.tim',skip=4)
        tmp <- tmp0[-c(1:4)]
        tmp <- gsub('.+ff|-fe.+','',tmp)
        data <- c()
        for(i in 1:length(tmp)){
            str <- unlist(strsplit(tmp[i],split=' '))
            str <- str[str!='']
            data <- rbind(data,as.numeric(as.character(str)))
        }
    }
}
t1 <- time_Jd2yr(OutTime$tB)
plot(t1,time_T2mT2(OutTime$tauE,OutTime$tB),xlab='tB [yr]',ylab='tauE-tB [s]',type='l')
plot(t1,time_T2mT2(OutTime$tB,OutTime$tS),xlab='tB [yr]',ylab='tB - tS [s]',type='l')
plot(t1,time_T2mT2(OutTime$BJDtcb,OutObs$JDtcb),xlab='tB [yr]',ylab='tS-JD[TCB] [s]',type='l')
plot(t1,time_T2mT2(OutObs$JDtcb,utc),xlab='tB [yr]',ylab='JD[TCB]-JD[UTC] [s]',type='l')

dev.off()

if(Par$binary & Par$Np>0){
fout <- paste0('../results/',Par$star,'_shapiro.pdf')
cat(fout,'\n')
pdf(fout,12,4)
par(mar=c(6,6,3,1),cex.lab=2,cex.axis=2,cex.main=2)
plot(phase[ind],OutTime$ShapiroTarget[ind]*1e6,xlab='Phase',ylab=expression('Shapiro Delay ['*mu*'s]'),type='l',main='PEXO')
dev.off()
}

if(FALSE){
xrange <-c(0,4)
t2 <- rowSums(utc-as.numeric(utc[1,]))
plot(t2,rowSums(OutObs$JDtt-utc)*DAYSEC,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
plot(t2,rowSums(OutObs$JDtcb-OutObs$JDtt)*DAYSEC,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
plot(t2,rowSums(tS-OutObs$JDtcb)*DAYSEC,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
plot(t2,rowSums(tS-utc)*DAYSEC,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')

if(Par$SBscaling){
plot(t1,OutTime$VacuumIS,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
plot(t1,OutTime$EinsteinIS,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
}
plot(t1,rowSums(OutTime$tB-OutTime$tS),xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')

if(Par$binary){
    plot(t1,(OutTime$RoemerTarget+OutTime$EinsteinTarget),xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
    plot(t1,OutTime$ShapiroTarget,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
    plot(t1,OutTime$AbeTarget,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
    plot(t1,rowSums(OutTime$tauE-OutTime$tB)*DAYSEC,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),type='l')
}

if(Par$SBscaling){
plot(t1[-1],diff(OutTime$VacuumIS)/diff(rowSums(tS)),xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),ylab=expression(dot(Delta)*t*' [second/day]'),type='l',xlim=xrange)
plot(t1[-1],diff(OutTime$EinsteinIS)/diff(rowSums(tS)),xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),ylab=expression(dot(Delta)*t*' [second/day]'),type='l',xlim=xrange)
}
plot(t1[-1],diff(Dt)/diff(rowSums(tS))*DAYSEC,xlab=expression('JY[UTC]-'*JY[0]*'[UTC]'),ylab=expression(dot(Delta)*t*' [second/day]'),type='l',xlim=xrange)
}


