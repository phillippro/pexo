source('general_function.R')
source('astrometry_function.R')
###see https://github.com/Starlink/starlink/blob/master/libraries/sla/refcoq.f
###for details of the test
hm <- 2380#height above sea level (meter)
tdk <- 280.15#ambient temperature at the observer (K)
pmb <- 1005#pressure at the observer (millibar)
rh <- 0.8#relative humidity at the observer (range 0-1)
wl <- 0.574#effective wavelength of the source (micrometer)
phi <- 50/180*pi#latitude of the observer (radian, astronomical)
tlr <- 0.0065#temperature lapse rate in the troposphere (K/meter)
eps <- 1e-8#precision required to terminate iteration (radian)
zreal.deg <- seq(5,85,by=1)
zreal <- zreal.deg/180*pi#observed zenith of the source
#zreal <- 10/180*pi#observed zenith of the source
Ntry <- 100
tol <- 1e-10
t0 <- proc.time()
zobs <- zreal

ref1 <- astro_refro(zreal,hm,tdk,pmb,rh,wl,phi,tlr,eps)
dur1 <- as.numeric(proc.time()-t0)[3]
cat('duration of refro:',dur1,'s\n')

t0 <- proc.time()
ind0 <- which(zreal.deg>=5 & zreal.deg<=85)
ind1 <- which(zreal.deg<5 | zreal.deg>85)
ref2 <- rep(NA,length(zreal.deg))
if(length(ind0)>0){
    refab <- astro_refco(hm,tdk,pmb,rh,wl,phi,tlr,eps)
    A <- refab$refa
    B <- refab$refb
    ref2[ind0] <- A*tan(zobs[ind0])+B*tan(zobs[ind0])^3
}
if(length(ind1)>0){
    ref2[ind1] <- astro_refro(zreal[ind1],hm,tdk,pmb,rh,wl,phi,tlr,eps)
}
dur2 <- as.numeric(proc.time()-t0)[3]
cat('duration of refco:',dur2,'s\n')

t0 <- proc.time()
ref3 <- rep(NA,length(zreal.deg))
if(length(ind0)>0){
    refab <- astro_refcoq(tdk,pmb,rh,wl)
    A <- refab$refa
    B <- refab$refb
    ref3[ind0] <- A*tan(zobs[ind0])+B*tan(zobs[ind0])^3
}
if(length(ind1)>0){
    ref3[ind1] <- astro_refro(zreal[ind1],hm,tdk,pmb,rh,wl,phi,tlr,eps)
}
dur3 <- as.numeric(proc.time()-t0)[3]
cat('duration of refcoq:',dur3,'s\n')

zd <- zobs/pi*180
out <- cbind(zd,ref1*206265,ref2*206265,ref3*206265)
colnames(out) <- c('ZD[deg]','refro[as]','refco[as]','refcoq[as]')
print(out)

fpdf <- '../results/compare_different_refraction_code.pdf'
cat(fpdf,'\n')
pdf(fpdf,6,4)
#par(mfrow=c(2,2))
par(mar=c(4,4,1,5))
durs <- c(dur1,dur2,dur3)
ns <- c('REFRO','REFCO','REFCOQ')
rad2am <- 180/pi*60
plot(zobs*180/pi,ref1*rad2am,xlab='z[deg]',ylab='R [arcmin]',type='l',ylim=range(ref1,ref2,ref3)*rad2am)
lines(zobs*180/pi,ref2*rad2am,col='red')
lines(zobs*180/pi,ref3*rad2am,col='blue')
legend('topright',inset=c(-0.1,0),legend=paste(ns,round(durs,2)),col=c('black','red','blue'),lty=1,xpd=NA,bty='n',horiz=TRUE)
dev.off()

