source('OrbitFunction.R')
library(magicaxis)
pdf('potential_mass.pdf',9,6)
size <- 1.2
par(mar=c(5,5,1,14),cex=size,cex.lab=size,cex.axis=size)
###constants
Me <- 0.000003003#Msol
Me.kg <- 5.972e24#kg
au2km <- 1.496e8
Msun.kg <- 1.98847e30#kg
Rsun.km <- 6.957e5#km
Re <- 6356#km
c <- 299792458#m/s
eta <- (4.74047e3/c)^2

##Pound & Rebka 1959
logphi.pr <- log10(4*pi^2*Me/(Re/au2km)*eta)
logmass.pr <- log10(Me.kg*1e3)

## precession of mercury
rmin.mer <- 0.307499#au
rmax.mer <-0.466 #au
a.mer <- 0.387#au
logphi.mer <- log10(4*pi^2/a.mer*eta)
mass.mer <- 3.3011e23#kg
logmass.mer <- log10(Msun.kg*1e3)

##light deflection and shapiro delay
logphi.len <- log10(4*pi^2/(Rsun.km/au2km)*eta)
logmass.len <- log10(Msun.kg*1e3)

##Hulse-Taylor pulsar
Mc <- 1.387
Mtot <- 2.828378
Mt <- Mtot-Mc
a.km <- 1950100# km
q.km <- 746600#km
logphi.ht <- log10(4*pi^2*Mtot/(q.km/au2km)*eta)
logmass.ht <- log10(Mtot*Msun.kg*1e3)

##gravitational redshift of Sirius B
dvr.gr <- 80.42#km/s
logphi.sb <- log10(dvr.gr*1e3/c)
logmass.sb <- log10(1.018*Msun.kg*1e3)

##S2 orbit around Sgr A*
z <- 200e3/c
logphi.s2 <- log10(z)
logmass.s2 <- log10(4.1e6*Msun.kg*1e3)

##LIGO detections
mlow <- 1.46+1.27
mup <- 50.6+34.3
logphi.ligo <- rep(log10(0.5),2)
logmass.ligo <- log10(c(mlow,mup)*Msun.kg*1e3)

##Event horizon
M87 <- 2400e12
logphi.eht <- log10(0.5)
logmass.eht <- log10(2400e12)+log10(Msun.kg*1e3)

##Fe alpha
#fabian00
logphi.fe <- log10(1/6)
logmass.fe <- 8+log10(Msun.kg*1e3)

##binary gravitational redshift
##alpha Cen
logphi.ac <- log10(0.25/c)
logmass.ac <- log10(Msun.kg*1e3)
Ms <- c(0.01,150)
vs <- c(0.01,1000)
logphi.binary <- log10(vs/c)
logmass.binary <- log10(Ms*Msun.kg*1e3)
#logphi.bi.low <- log10(0.1/c)


####milky way dark matter
Mdisk <- 6.8e10
Mnuc <- 0.17e10
Mbul <- 0.5e10

##ref. Lelli 2017; Milky-way like galaxy
gcrit <- gup0 <- 1e-10#m/s^2
pc2km <- 3.086e13#km
Rcrit <- 5#kpc
v2 <- (gcrit*pc2km)*Rcrit*1e3*1e3#m^2/s^2
vcrit <- sqrt(v2)#m/s
logphi.rot <- log10((vcrit/c)^2)
logmass.rot <- log10((vcrit/1e3/4.74047)^2*Rcrit*1e3*206265/(4*pi^2))+log10(Msun.kg*1e3)


####upper boundary at g=glow
#Rcs <- seq(0.1,1000,by=0.1)
#v2 <- (gcrit*pc2km)*Rcs*1e3*1e3#m^2/s^2
Ms <- c(1e-3,1e23)
gup <- gup0/4.74047/1e3*365.25*3600*24#m/s^2 to au/yr^2
r <- sqrt(4*pi^2*Ms/gup)#au
v2 <- 4*pi^2*Ms/r
logphi.rots1 <- log10(v2*eta)
#logmass.rots1 <- log10((vcrit/1e3/4.74047)^2*Rcs*1e3*206265/(4*pi^2))+log10(Msun.kg*1e3)
logmass.rots1 <- log10(Ms)+log10(Msun.kg*1e3)

####lower boundary at g=gup
glow0 <- 1e-12
#Rcs <- seq(1e-4,1e4,by=0.1)#from 0.1 Msol to 10^14 Msol
#v2 <- (gcrit*pc2km)*Rcs*1e3*1e3#m^2/s^2
glow <- glow0/4.74047/1e3*365.25*3600*24#m/s^2 to au/yr^2
r <- sqrt(4*pi^2*Ms/glow)#au
v2 <- 4*pi^2*Ms/r
logphi.rots2 <- log10(v2*eta)
logmass.rots2 <- log10(Ms)+log10(Msun.kg*1e3)


##wide binary
a.proxima <- 8700#au
q.proxima <- 8700*(1-0.50)#au
r.proxima <- 13e3#au
Mab <- 2.0429#Msol
v2 <- 4*pi^2*Mab/r.proxima
logphi.wb <- log10(v2*eta)
logmass.wb <- log10(Mab*Msun.kg*1e3)

##wide binaries
ms <- seq(0.1,1)
r.low <- 3e3#au
r.up <- 206265#au
a.wb.up <- 4*pi^2/r.low^2
##phi determined by v contour
logphi.wbv.up <- log10(4*pi^2/r.low)
logphi.wbv.low <- log10(4*pi^2/r.up)
##phi determined by a contour
a.wb.low <- 4*pi^2/r.up^2
logmass

##micro-lensing
#https://iopscience.iop.org/article/10.1086/178096/pdf
Ms.len <- c(1e-7,1e14)
D <- c(1e-3,1e0,1e6)
logphi.close <- log10(4.03*sqrt(Ms*D[1]/2)/(2*D[1])/206265e3)
logphi.av <- log10(4.03*sqrt(Ms*D[2]/2)/(2*D[2])/206265e3)
logphi.far <- log10(4.03*sqrt(Ms*D[3]/2)/(2*D[3])/206265e3)
logmass.close <- logmass.av <- logmass.far <- log10(Ms.len)+log10(Msun.kg*1e3)

##all combined
logphis <- c(logphi.pr,logphi.mer,logphi.len,logphi.ht,logphi.sb,logphi.s2,logphi.ligo,logphi.eht,logphi.fe,logphi.ac,logphi.wb)
logmass <- c(logmass.pr,logmass.mer,logmass.len,logmass.ht,logmass.sb,logmass.s2,logmass.ligo,logmass.eht,logmass.fe,logmass.ac,logmass.wb)
ns <- c('pr','mer','len','ht','sb','s2','ligo.low','ligo.up','eht','fe','ac','wb')
names(logmass) <- names(logphis) <- ns

plot(logmass,logphis,xlab=expression(log[10]*'(mass/gr)'),ylab=expression(log[10]*'('*phi/c^2*')'),pch=20,cex=2.5,ylim=c(min(logphis)-0.5,max(logphis)+1),xaxt='n',yaxt='n')
xoff <- 2
yoff <- 1.8
magaxis(side=1:2)
points(logmass['ac'],logphis['ac'],col='red',pch=20,cex=2.5)
#text(logmass,logphis,labels=ns,pos=1,xpd=TRUE)

xmax <- max(logmass)+xoff
####label each experiment
off <- -1.8
lines(c(logmass.ligo[1],xmax),c(logphi.ligo[1],off),xpd=NA,lwd=2,col='grey')
lines(c(logmass.ligo[2],xmax),c(logphi.ligo[2],off),xpd=NA,lwd=2,col='grey')
text(xmax,off,label='GW Detections by LIGO',pos=4,xpd=TRUE)

lines(c(logmass.eht,xmax),c(logphi.eht,logphi.eht),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.eht,label='EHT imaging',pos=4,xpd=TRUE)

off <- -0.3
lines(c(logmass.fe,xmax),c(logphi.fe,logphi.fe+off),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.fe+off,label=expression('Fe K-'*alpha*' line'),pos=4,xpd=TRUE)

off <- 0
lines(c(logmass.s2,xmax),c(logphi.s2,logphi.s2+off),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.s2+off,label=expression('S2 orbit around Sgr'*A^{'*'}),pos=4,xpd=TRUE)

off <- -0.2
lines(c(logmass.sb,xmax),c(logphi.sb,logphi.sb+off),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.sb+off,label='Self grav. redshift, Sirius B',pos=4,xpd=TRUE)

off <- 0.5
lines(c(logmass.ht,xmax),c(logphi.ht,logphi.ht+off),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.ht+off,label='Hulse-Taylor pulsar',pos=4,xpd=TRUE)

lines(c(logmass.len,xmax),c(logphi.len,logphi.len),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.len,label='Light deflection & Shapiro delay',pos=4,xpd=TRUE)
#text(xmax,logphi.len-0.3,label='Shapiro Delay',pos=4,xpd=TRUE)
#lines(c(xmax,xmax),c(logphi.len+0,logphi.len-0.5),xpd=NA,lwd=2,col='grey')

lines(c(logmass.rot,xmax),c(logphi.rot,logphi.rot),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.rot,label='Milky way rotation curve',pos=4,xpd=TRUE)

lines(c(logmass.mer,xmax),c(logphi.mer,logphi.mer),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.mer,label='Mercury precession',pos=4,xpd=TRUE)

lines(c(logmass.ac,xmax),c(logphi.ac,logphi.ac-0.6),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.ac-0.6,label=expression(alpha*' Centauri gravitational redshift'),pos=4,xpd=TRUE)
#
polygon(logmass.binary[c(1,2,2,1)],logphi.binary[c(1,1,2,2)],col=tcol('red',50),border=FALSE)
lines(c(logmass.ac+1,xmax),c(logphi.ac,logphi.ac)+0.5,xpd=NA,lwd=2,col='grey')
text(xmax,logphi.ac+0.5,label=expression('Binary gravitational redshift'),pos=4,xpd=TRUE)

lines(c(logmass.pr,xmax),c(logphi.pr,logphi.pr-1.3),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.pr-1.3,label='Pound-Rebka experiment',pos=4,xpd=TRUE)

points(logmass.rot,logphi.rot,pch=20,cex=2.5,col='darkgrey')
#lines(logmass.rots1,logphi.rots1,col='darkgrey')
#lines(logmass.rots1,logphi.rots1,col='darkgrey')
polygon(c(logmass.rots1,rev(logmass.rots2)),c(logphi.rots1,rev(logphi.rots2)),col=tcol('grey',50),border=FALSE)

##wide binary
points(logmass.wb,logphi.wb,pch=20,cex=2.5,col='darkgrey')
lines(c(logmass.wb,xmax),c(logphi.wb,logphi.wb-0.5),xpd=NA,lwd=2,col='grey')
text(xmax,logphi.wb-0.5,label=expression('Proxima Centauri acceleration'),pos=4,xpd=TRUE)

##rotation or wide binary
lines(c(logmass.wb+2,xmax),c(logphi.wb,logphi.wb)+0.3,xpd=NA,lwd=2,col='grey')
text(xmax,logphi.wb+0.3,label=expression('Galactic rotation or wide binary'),pos=4,xpd=TRUE)

##Lensing at Einstein radius
#lines(logmass.close,logphi.close,lty=2)
#lines(logmass.av,logphi.av,lty=2)
#lines(logmass.far,logphi.far,lty=2)
#polygon(logmass.close[c(1,2,2,1)],c(logphi.close[1],min(logphi.close[2],0),min(logphi.far[2],0),logphi.far[1]),col=tcol('blue',50),border=FALSE)
polygon(logmass.close[c(1,2,2,1)],c(logphi.close,rev(logphi.far)),col=tcol('blue',60),border=FALSE)
lines(c(39,xmax),c(-2.5,-2.5),xpd=NA,lwd=2,col='grey')
text(xmax,-2.5,label=expression('Lensing from 1 pc to 1 Gpc'),pos=4,xpd=TRUE)
#text(xmax,logphi.close+0.3,label=expression('Lensing at D=1 kpc'),pos=4,xpd=TRUE)
#text(xmax,logphi.close+0.3,label=expression('Lensing at D=1 Mpc'),pos=4,xpd=TRUE)

### black hole horizon
polygon(x=c(min(logmass)-xoff,max(logmass)+xoff,max(logmass)+xoff,min(logmass)-xoff),y=c(0,0,yoff,yoff),col='black')
text(0.5*(min(logmass)+max(logmass)),0.5,label='Black Holes',col='white')
dev.off()
