source('read_TEMPO2output.R',local=TRUE)
if(UNITS=='TCB'){
    F0 <- 1
}else{
#    F0 <- 1.00000001550505
    F0 <- 1
}
#F0 <- 1.00000001550505
f0 <- 5.45e8#HZ
DJM0 <- 2400000.5
CMPS <- 299792458.0
pexo <- read.table('../data/10700freqPEXO.dat',header=TRUE)
t2 <- read.table('../data/10700freqT2.dat',header=TRUE)

###TEMPO2 prediction/emulation code
if(FALSE){
    if(UNITS=='TCB'){
        system("source ~/Documents/tempo2/emulation_tcb.sh")
    }else{
        system("source ~/Documents/tempo2/emulation_tdb.sh")
    }
}
te <- read.table('../data/tempo2_emulation.dat',header=TRUE)
fname <- '../results/compare_frequency.pdf'
#cat('Frequency test results:\n',fname,'\n')
cat('Frequency test results:\n')
cat('UNITS:',UNITS,'\n')
BJDtcb <- t2[1,'BJDtcb']+DJM0
#cat('general2 output z=',(545049587682034.5-545000000e6)/545000000e6,'\n')
mbjd.tdb <- (bjd.tdb[,1]-DJM0)+bjd.tdb[,2]
#cat('mbjd.tdb=',mbjd.tdb,'\n')
mbjd.tcb <- (bjd.tcb[,1]-DJM0)+bjd.tcb[,2]
mjd.tt <- (tb$tt[,1]-DJM0)+tb$tt[,2]
mjd.tcb <- (tb$tcb[,1]-DJM0)+tb$tcb[,2]
mjd.tdb <- (tb$tdb[,1]-DJM0)+tb$tdb[,2]
mjd.utc <- (utc[,1]-DJM0)+utc[,2]
mbjd.tdb <- time_Jd2mjd(bjd.tdb)
mbjd.tcb <- time_Jd2mjd(bjd.tcb)
#cat('mbjd.tcb=',mbjd.tcb,'\n')
f0 <- 545000000e6
if(UNITS=='TCB'){
#    fSSB <- c(545049575262787.4375,545049694157905.5)
#    fSSB <- c(545049575262787.4375,545049694158471.1875)#10700w.tim
    fSSB <- c(545049575262787.4375,544998913503703.875,544952061593343.0625,545014131689073.25)#10700w1.tim
#    cat('general2 output z=',(545049579231701.1875-545000000e6)/545000000e6,'\n')
}else{
    fSSB <- 545049583713120.9375
}
zg <- (fSSB-f0)/f0
cat('general2 output z=',zg,'\n')
#mjd <- mjd.tt
#mjd <- mjd.tcb
mjd <- mjd.tdb
#mjd <- mbjd.tcb
#mjd2 <- cbind(tb$tt[,1]-DJM0,tb$tt[,2])
#mjd2 <- cbind(tb$tcb[,1]-DJM0,tb$tcb[,2])
#mjd2 <- cbind(tb$tdb[,1]-DJM0,tb$tdb[,2])
mjd2 <- cbind(bjd.tdb[,1]-DJM0,bjd.tdb[,2])
#mjd2 <- cbind(bjd.tcb[,1]-DJM0,bjd.tcb[,2])
ind <- sapply(1:length(mbjd.tcb), function(i) which(mjd[i]<te[,'mjd'])[1]-1)
#ind <- sapply(1:length(mjd), function(i) which.min(abs(mjd[i]-te[,'mjd'])))
#    dt <- (mbjd.tcb-te[ind,'mjd'])*1440
dt <- ((mjd2[,1]-te[ind,'mjd'])+mjd2[,2])*1440
dF <- rowSums(sapply(2:12,function(i) (i-1)*te[ind,paste0('coeff',i)]*dt^(i-2)))/60
ze <- F0/(F0+dF)-1#
ze1 <- -dF
ze2 <- -dF+dF^2
cat('polyco emulation output ze=',ze,'\n')
cat('polyco emulation output ze1=',ze1,'\n')
cat('polyco emulation output ze2=',ze2,'\n')
cat('PEXO output =',rv$zcomb[,'obs'],'\n')
if(FALSE){
    z <- rv$zcomb
    zsep <- z[,'gvo']-z[,'ro']-z[,'so']
    cat('PEXOsep output: gvo-ro-so =',zsep,'\n')
    cat('RV PEXOsep-general2 =',(zsep-zg)*CMPS,'m/s\n')
    cat('RV PEXO-general2 =',(rv$zcomb[,'obs']-zg)*CMPS,'m/s\n')
}

#####plot
if(nrow(te)==nrow(t2)){
pdf(fname,8,8)
par(mfrow=c(2,2))
zp <- pexo[,'zobs']
vp <- zp*CMPS
#vp2 <- te[,'doppler']*1e-4*CMPS
ve <- ze*CMPS
zt <- (t2[,'fSSB']-t2[,'fOBS']*1e6)/(t2[,'fOBS']*1e6)
vt <- zt*CMPS
plot(BJDtcb,vp)
plot(BJDtcb,vt)
plot(BJDtcb,vp-vt)
plot(te[,'mjd'],vp2)
fte <- approxfun(te[,'mjd'],vp2)
#vp3 <- fte(t2[,'BJDtcb'])

#plot(vp,vt)
dev.off()
}
