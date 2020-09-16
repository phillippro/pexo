source('OrbitFunction.R')
rvA <- OutRv$HD239960$rv$RvST
rvB <- OutRv$HD239960$rv$RvSB-OutRv$HD239960$rv$RvBT*Par$mT/exp(ParOpt['logmC1'])
fout <- '../results/kruger60.pdf'
cat(fout,'\n')
pdf(fout,8,4)
par(mfrow=c(1,2))
jd1 <- rowSums(OutObs$HD239960$rel$JDutc)
jd2 <- rowSums(OutObs$HD239960$rv$JDutc)
rv1 <- OutRv$HD239960$rel$RvSB
#dt1 <- (jd1-Par$epoch)/DJY
epoch <- 2453000
dt1 <- (jd1-epoch)/DJY
tmp <- lm(rv1~dt1)
y1 <- tmp$coefficients[1]+tmp$coefficients[2]*dt1
plot(jd1,rv1,xlab='JD [day]',ylab='systematic rv [m/s]',xlim=range(jd1,Par$epoch),type='l',ylim=range(rv1,y1),col='white')
lines(jd1,y1)
abline(v=epoch,col='green')
abline(v=Par$epoch,col='blue')
cat('rv(t0)=',tmp$coefficients[1],'m/s\n')
rv1 <- tmp$coefficients[1]#at a given reference epoch
legend('top',xpd=NA,inset=c(0,-0.2),legend=c('Gaia epoch','KECK epoch'),col=c('blue','green'),lty=1,bty='n',horiz=TRUE)

plot(jd2,rvA,xlab='JD [day]',ylab='Binary RV [m/s]',ylim=range(c(rvA,rvB)))
points(jd2,rvB,col='red')
rv0 <- (ParOpt['rvOff']+Par$rv)*1e3
cat('rv0=',rv0,'m/s\n')
abline(h=rv0,col='blue')
abline(h=rv1,col='green')
legend('top',xpd=NA,inset=c(0,-0.2),legend=c('A','B'),col=c('black','red'),pch=1,bty='n',horiz=TRUE)

#plot(c(-1,1),c(-1,1),col='white',axes=FALSE,xlab='',ylab='')
x <- seq(-0.8,0.8,length.out=4)
y <- seq(1,-1,by = -0.2)

ra0 <- Par$ra
dra <- ParOpt['raOff']*DMAS2R
ra <- ra0+dra

dec0 <- Par$dec
ddec <- ParOpt['decOff']*DMAS2R
dec <- ddec+dec0

plx0 <- Par$plx
dplx <- ParOpt['plxOff']
plx <- dplx+plx0

pmra0 <- Par$pmra
dpmra <- ParOpt['pmraOff']
pmra <- pmra0+dpmra

pmdec0 <- Par$pmdec
dpmdec <- ParOpt['pmdecOff']
pmdec <- pmdec0+dpmdec

rv0 <- Par$rv
drv <- ParOpt['rvOff']
rv <- rv0+drv

ns <- c('ra','dec','plx','pmra','pmdec','rv')
val1 <- c(ra0*180/pi,dec0*180/pi,plx0,pmra0,pmdec0,rv0)
names(val1) <- paste0(ns,'.gaia')
unit1 <- c('deg','deg','mas','mas/yr','mas/yr','km/s')

val2 <- c(dra/DMAS2R,ddec/DMAS2R,dplx,dpmra,dpmdec,drv)
names(val2) <- paste0(ns,'.offset')
unit2 <- c('mas','mas','mas','mas/yr','mas/yr','km/s')

val3 <- c(ra*180/pi,dec*180/pi,plx,pmra,pmdec,rv)
names(val3) <- paste0(ns,'.new')
unit3 <- c('deg','deg','mas','mas/yr','mas/yr','km/s')

for(k in 1:length(ns)){
    cat(names(val1)[k],'=',val1[k],unit1[k],';',names(val2)[k],'=',val2[k],unit2[k],';',names(val3)[k],'=',val3[k],unit3[k],'\n')
}

###SB in equatorial coordinates
d <- as.numeric(1/plx)#kpc
rSB <- as.numeric(bl2xyz(dec,ra))*d*1e3#pc
cat('Heliocentric position of Kruger 60 barycenter in equatorial coordiante system in units of pc:',rSB,'\n')
cat('Heliocentric position of Kruger 60 barycenter in Galactic coordiante system in units of pc:',e2g.vel(rSB),'\n')

vSB <- OutTime$HD239960$rv$vSB[1,]*auyr2kms
cat('Heliocentric velocity of Kruger 60 barycenter in equatorial coordiante system in units of km/s:',vSB,'\n')
cat('Heliocentric velocity of Kruger 60 barycenter in Galactic coordiante system (UVW) in units of km/s:',e2g.vel(vSB),'\n')
dev.off()
