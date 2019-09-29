dir <- '../results/'
fname <- paste0(dir,star,'_df',df,'_tbase',tbase,'_dt',dt,'order',Norder,'type',type,'.pdf')
cat('output pdf:\n',fname,'\n')
pdf(fname,8,8)
par(mfrow=c(2,2),mar=c(5,5,1,1))
dv <- (zd-zg)*CMPS*1e3
ind.out <- which(abs(dv)>0.2)
JDutc <-rowSums(utc1)
tout <- JDutc[ind.out]
plot(JDutc,dv,xlab='JD[UTC]',ylab=expression(Delta*v[r]^{TEMPO2}-Delta*v[r]^{PEXO}*' [mm/s]'),type='l')

plot(JDutc,tb$tt[index,2]-utc[index,2],xlab='JD[UTC]',ylab=expression(Delta*v[P]-Delta*v[T]),type='l')

ii <- 1:nrow(utc1)
jj <- 1:nrow(utc)
plot(JDutc[ii],bjd.tcbs[index[ii],2]-bjd.tcb[index[ii],2],xlab='JD[UTC]',ylab='BJD[TCB] (TEMPO2)')
plot(JDutc[ii],jd.tts[index[ii],2]-tb$tt[index[ii],2],xlab='JD[UTC]',ylab='JD[TT] (TEMPO2)')
dt1 <- jd.tts[index[ii],2]-tb$tt[index[ii],2]
dt2 <- jd.tts[index[ii]+1,2]-tb$tt[index[ii]+1,2]
dt3 <- jd.tts[index[ii]-1,2]-tb$tt[index[ii]-1,2]
ddt1 <- (jd.tts[index[ii]+1,2]-jd.tts[index[ii]-1,2])-(tb$tt[index[ii]+1,2]-tb$tt[index[ii]-1,2])
ddt2 <- (bjd.tcbs[index[ii]+1,2]-bjd.tcbs[index[ii]-1,2])-(bjd.tcb[index[ii]+1,2]-bjd.tcb[index[ii]-1,2])
dbjd.tempo <- (bjd.tcbs[index[ii]+1,2]-bjd.tcbs[index[ii]-1,2])
dbjd.pexo <- (bjd.tcb[index[ii]+1,2]-bjd.tcb[index[ii]-1,2])
dbjd1 <- bjd.tcbs[index[ii]-1,2]-bjd.tcb[index[ii]-1,2]
dbjd2 <- bjd.tcbs[index[ii],2]-bjd.tcb[index[ii],2]
dbjd3 <- bjd.tcbs[index[ii]+1,2]-bjd.tcb[index[ii]+1,2]
plot(JDutc[ii],dt1,xlab='JD[UTC]',ylab='dt123',ylim=range(dt1,dt2,dt3))
points(JDutc[ii],dt2,col='blue')
points(JDutc[ii],dt3,col='red')
plot(JDutc[ii],ddt1,xlab='JD[UTC]',ylab='ddt1')
plot(JDutc[ii],ddt2,xlab='JD[UTC]',ylab='ddt2')
#abline(v=tout,col='red')
plot(JDutc[ii],dbjd.tempo,xlab='JD[UTC]')
#abline(v=tout,col='red')
plot(JDutc[ii],dbjd.pexo,xlab='JD[UTC]')
#abline(v=tout,col='red')
plot(JDutc[ii],dbjd1,xlab='JD[UTC]')
#abline(v=tout,col='red')
plot(JDutc[ii],dbjd2,xlab='JD[UTC]')
#abline(v=tout,col='red')
plot(JDutc[ii],dbjd3,xlab='JD[UTC]')
#abline(v=tout,col='red')

plot(rowSums(utc)[jj],bjd.tcbs[jj,2]-bjd.tcb[jj,2],xlab='JD[UTC]',ylab='dBJD[TCB] (TEMPO2)')
plot(rowSums(utc)[jj],jd.tts[jj,2]-tb$tt[jj,2],xlab='JD[UTC]',ylab='dJD[TT] (TEMPO2)')
points(JDutc[ii],jd.tts[index[ii],2]-tb$tt[index[ii],2],col='blue')


plot(JDutc,tb$tt[index,2]-jd.tts[index,2],xlab='JD[UTC]',ylab='dJD[TT]')
plot(JDutc,tb$tt[index,2],xlab='JD[UTC]',ylab='JD[TT] (PEXO)')
plot(JDutc,jd.tts[index,2],xlab='JD[UTC]',ylab='JD[TT] (TEMPO2)')
plot(JDutc,tb$tai[index,2],xlab='JD[UTC]',ylab='JD[TAI] (PEXO)')

plot(JDutc,bjd.tcb[index,2]-bjd.tcbs[index,2],xlab='JD[UTC]',ylab='dBJD')
plot(JDutc,bjd.tcb[index,2],xlab='JD[UTC]',ylab='BJD[TCB] (PEXO)')
plot(JDutc,bjd.tcbs[index,2],xlab='JD[UTC]',ylab='BJD[TCB] (TEMPO2)')


###shapiro delays
dbjd <- ((bjd.tcbs[index+1,2]-bjd.tcbs[index-1,2])+(bjd.tcbs[index+1,1]-bjd.tcbs[index-1,1]))*DAYSEC
zs <- (data[index+1,7]-data[index-1,7])/dbjd
zsJ <- (data[index+1,8]-data[index-1,8])/dbjd
zsS <- (data[index+1,9]-data[index-1,9])/dbjd
zsV <- (data[index+1,10]-data[index-1,10])/dbjd
zsU <- (data[index+1,11]-data[index-1,11])/dbjd
zsN <- (data[index+1,12]-data[index-1,12])/dbjd
plot(JDutc,data[index,7],xlab='JD[UTC]',ylab='ShapiroSun')
plot(JDutc,zs*CMPS*1e3,xlab='JD[UTC]',ylab='ShapiroSun rv [mm/s]')
#abline(v=tout,col='red')
plot(JDutc,data[index,8],xlab='JD[UTC]',ylab='ShapiroJupiter')
plot(JDutc,zsJ*CMPS*1e3,xlab='JD[UTC]',ylab='ShapiroJ rv [mm/s]')
#abline(v=tout,col='red')
plot(JDutc,data[index,9],xlab='JD[UTC]',ylab='ShapiroSaturn')
plot(JDutc,zsS*CMPS*1e3,xlab='JD[UTC]',ylab='ShapiroS rv [mm/s]')
#abline(v=tout,col='red')
plot(JDutc,data[index,10],xlab='JD[UTC]',ylab='ShapiroVenus')
plot(JDutc,zsV*CMPS*1e3,xlab='JD[UTC]',ylab='ShapiroV rv [mm/s]')
#abline(v=tout,col='red')
plot(JDutc,data[index,11],xlab='JD[UTC]',ylab='ShapiroUranus')
plot(JDutc,zsU*CMPS*1e3,xlab='JD[UTC]',ylab='ShapiroU rv [mm/s]')
#abline(v=tout,col='red')
plot(JDutc,data[index,12],xlab='JD[UTC]',ylab='ShapiroNeptune')
plot(JDutc,zsN*CMPS*1e3,xlab='JD[UTC]',ylab='ShapiroN rv [mm/s]')
#abline(v=tout,col='red')

#plot(utc1,data[index,'tropo'],xlab='JD[UTC]',ylab='tropo',type='l')
dev.off()
