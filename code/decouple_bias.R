library(magicaxis)
library(fields)
cols <- rainbow(100,start=0.8)#
yy <- ms <- exp(seq(log(1e-4),log(10),by=0.1))
Ms <- ms+1#total mass
xx <- as <- exp(seq(log(1e-3),log(1e3),by=0.1))
ds <- c(1, 10,100,1000)
plxs <- 1000/ds
vtot <- 5e4#m/s
pdf('decouple_bias.pdf',16,16)
size <- 1.2
par(mfrow=c(2,2))
for(k in 1:length(ds)){
par(mar=c(6,6,6,8),cex=size,cex.lab=size,cex.axis=size)
    plx  <- plxs[k]
    zz <- array(NA,dim=c(length(xx),length(yy)))
    for(j in 1:length(as)){
        Ps <- sqrt(as[j]^3/Ms)#yr
        aa <- as[j]*ms/Ms#stellar reflex orbital semi-major axis
        vs <- 2*pi*aa/Ps#au/yr; stellar reflext velocity
        dpm <- vs*plx#mas/yr
        zz[j,] <- vtot*dpm/1e3/206265
    }
    zlim <- range(zz)
    image(log10(xx),log10(yy),zz,xlab='Semi-major axis [au]',ylab='Mass [Msun]',xaxt='n',yaxt='n',col=cols,main=paste0(ds[k],' pc'))
    magaxis(side=1:2,unlog=TRUE,tcl=-0.8)
    image.plot(zz,col=cols,legend.only=TRUE,zlim = zlim,legend.mar=7.1)
}
dev.off()
####
pdf('decouple_simple.pdf',6,6)
size <- 1.2
par(mfrow=c(1,1),mar=c(5,5,3,3),cex=size,cex.lab=size,cex.axis=size)
for(k in 1:length(ds)){
    plx <- 1000/ds[k]#mas
    vs.sqa <- 2*pi*ms/sqrt(Ms)
    tmp <- 1e-2/vtot*1e3*206265/plx/vs.sqa
    as <- 1/tmp^2#1 cm/s/yr threshold
    tmp2<- 2e-3/vtot*1e3*206265/plx/vs.sqa
    as2 <- 1/tmp2^2#2 mm/s/yr threshold
    as3 <- 1e-2/vtot*1e3*206265/plx*Ms/ms#1 cm/s threshold
    as4 <- 2e-3/vtot*1e3*206265/plx*Ms/ms#1 cm/s threshold
    if(k==1){
        plot(as,ms,xlab='Semi-major axis [au]',ylab='Mass [Msun]',xlim=range(xx),ylim=c(min(yy),1),type='l',log='xy',xaxt='n',yaxt='n')
#        magaxis(side=1:2,unlog=TRUE,tcl=-0.8)
        magaxis(side=1:4,unlog=TRUE)
    }else{
        lines(as,ms)
    }
    lines(as2,ms,lty=2)
    lines(as3,ms,col='grey')
    lines(as4,ms,lty=2,col='grey')
}
dev.off()

pdf('barycorr_MP.pdf',6,6)
size <- 1.2
Mj2s <- 9.543e-4#Jupiter mass in unit of solar mass
Mj2e <- 317.8#Jupiter mass in unit of Earth mass
par(mfrow=c(1,1),mar=c(5,5,3,3),cex=size,cex.lab=size,cex.axis=size)
for(k in 1:length(ds)){
    plx <- 1000/ds[k]#mas
    vs.sqa <- 2*pi*ms/sqrt(Ms)
    tmp <- 1e-2/vtot*1e3*206265/plx/vs.sqa
    as <- 1/tmp^2#1 cm/s threshold
    Ps <- sqrt(as^3/Ms)*365.25#day
    ms2 <- ms/Mj2s#Mjup
#    tmp2<- 2e-3/vtot*1e3*206265/plx/vs.sqa
#    as2 <- 1/tmp2^2#1 cm/s threshold
    if(k==1){
        plot(Ps,ms2,xlab='Orbital Period [day]',ylab='Mass [Msun]',type='l',log='xy',xaxt='n',yaxt='n')#ylim=min(yy/Mj2s)
#        magaxis(side=1:2,unlog=TRUE,tcl=-0.8)
        magaxis(side=1:4,unlog=TRUE)
    }else{
        lines(Ps,ms2)
    }
}
dev.off()
